"""Declared higher-link follower-force survey with independently checked fields."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from scipy.stats import qmc
from follower_link import matrices


def field(theta,p):
    n=len(theta);m=n-1
    stiffness=np.r_[1.,p[1:1+m]];drag=np.r_[1.,p[1+m:1+2*m]];length=np.r_[1.,p[1+2*m:1+3*m]]
    resistance,torque=matrices(theta,p[0],stiffness,drag,length)
    difference=np.eye(n)-np.eye(n,k=-1)
    elastic=difference.T@np.diag(stiffness)@difference
    rest=np.r_[0.,np.cumsum(p[1+3*m:])]
    return np.linalg.solve(resistance,torque+elastic@rest)


def write(path,record):
    path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    n=args.links;m=n-1;assert n>=2
    kernel=Path(__file__).with_name('follower_chain.cl')
    device=cl.get_platforms()[0].get_devices()[0];context=cl.Context([device]);queue=cl.CommandQueue(context)
    program=cl.Program(context,kernel.read_text()).build(options=[f'-D N={n}'])
    advance=cl.Kernel(program,'advance');evaluate=cl.Kernel(program,'evaluate')
    u=qmc.Sobol(1+4*m,scramble=True,seed=17400+n).random_base2(args.power)
    lo=np.r_[.6,np.full(m,.3),np.full(m,.5),np.full(m,.7)]
    hi=np.r_[100.,np.full(m,3.),np.full(m,2.),np.full(m,1.4)]
    params=np.column_stack((np.exp(np.log(lo)+u[:,:1+3*m]*np.log(hi/lo)),.4*u[:,1+3*m:]-.2))
    # First half has zero rest curvature; second half has small signed rest turns.
    params[:len(u)//2,1+3*m:]=0.
    initial=np.zeros((len(u),2*n));initial[:,1:n]=np.cumsum(params[:,1+3*m:],axis=1)
    initial[:,:n]+=np.linspace(.1,-.03,n);initial[:,n]=1.
    states=initial.copy();sums=np.zeros(len(u));status=np.zeros(len(u),np.int32);mf=cl.mem_flags
    sb=cl.Buffer(context,mf.READ_WRITE|mf.COPY_HOST_PTR,hostbuf=states)
    pb=cl.Buffer(context,mf.READ_ONLY|mf.COPY_HOST_PTR,hostbuf=params)
    lb=cl.Buffer(context,mf.READ_WRITE|mf.COPY_HOST_PTR,hostbuf=sums)
    fb=cl.Buffer(context,mf.READ_WRITE|mf.COPY_HOST_PTR,hostbuf=status)
    vb=cl.Buffer(context,mf.WRITE_ONLY,states.nbytes);values=np.empty_like(states)
    evaluate(queue,(len(u),),None,sb,pb,vb);cl.enqueue_copy(queue,values,vb).wait();worst=0.
    for i in range(len(u)):
        f=field(states[i,:n],params[i]);jv=field(states[i,:n]+1e-25j*states[i,n:],params[i]).imag/1e-25
        worst=max(worst,float(np.max(abs(values[i]-np.r_[f,jv]))/(1+np.max(abs(np.r_[f,jv])))))
    assert worst<1e-10
    np.savez_compressed(out/'declared-inputs.npz',parameters=params,initial=initial)
    write(out/'protocol.json',dict(links=n,count=len(u),seed=17400+n,step=args.step,warm=args.warm,duration=args.duration,
        logarithmic_bounds=[lo.tolist(),hi.tolist()],rest_turn_bounds=[-.2,.2],straight_count=len(u)//2,
        parameters='force, k[2:N], drag[2:N], length[2:N], rest_turn[2:N]',
        method='FP64 RK4 state plus analytic tangent, .1 tau renormalization; 4 observed windows',device=device.name,
        field_relative_error=worst,kernel_sha256=hashlib.sha256(kernel.read_bytes()).hexdigest(),
        sources=['De Canio2017 Sec.4','Ling, Guo, Kanso2018 Sec.7, DOI10.1098/rsif.2018.0594'],
        limitations='Established chain construction with declared parameter variations. Strict planar constraint; no contact or nonlocal drag. No novelty claim.'))
    for name in ['follower_chain.cl','follower_link.py','search_follower_chain.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    steps=int(round(25/args.step));qr=int(round(.1/args.step));assert steps%qr==0
    assert args.warm%25==0 and args.duration%100==0
    started=time.perf_counter();elapsed=0.;previous=np.zeros(len(u));windows=[]
    while elapsed<args.warm+args.duration-1e-9:
        advance(queue,(len(u),),None,sb,pb,lb,fb,np.float64(args.step),np.int32(steps),np.int32(qr)).wait();elapsed+=25.
        if elapsed==args.warm:sums[:]=0.;cl.enqueue_copy(queue,lb,sums).wait()
        if elapsed>args.warm and abs((elapsed-args.warm)%(args.duration/4))<1e-9:
            cl.enqueue_copy(queue,states,sb);cl.enqueue_copy(queue,sums,lb);cl.enqueue_copy(queue,status,fb).wait()
            windows.append((sums-previous)/(args.duration/4));previous=sums.copy()
            np.savez_compressed(out/f'checkpoint-{int(elapsed):06}.npz',states=states,sums=sums,status=status,windows=windows,elapsed=elapsed)
        if elapsed%100==0:print('tau',elapsed,'wall_seconds',round(time.perf_counter()-started,1),flush=True)
    rates=sums/args.duration;windows=np.array(windows).T
    keep=(status==0)&np.isfinite(rates)&(windows.min(axis=1)>.01)
    np.savez_compressed(out/'complete.npz',parameters=params,initial=initial,states=states,rate=rates,windows=windows,status=status,retained=keep)
    write(out/'analysis.json',dict(completed=int((status==0).sum()),failed=int((status!=0).sum()),retained_ids=np.flatnonzero(keep).tolist(),
        retained_rates=rates[keep].tolist(),wall_seconds=time.perf_counter()-started,
        retention='All four windows > .01 per tau and no numerical failure. Needs independent adaptive checks and contact audit.'))
    print('retained',int(keep.sum()),'of',len(u),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--links',type=int,default=4);p.add_argument('--power',type=int,default=11)
    p.add_argument('--step',type=float,default=.002);p.add_argument('--warm',type=float,default=200.)
    p.add_argument('--duration',type=float,default=1000.);main(p.parse_args())

"""GPU pilot of a physically declared, curved three-link extension.

No model novelty is asserted: follower-force filaments and intrinsic curvature
are established. This checks geometry only after a reproducible numerical gate.
"""
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
    stiffness=np.array([1.,p[1],p[2]])
    drag=np.array([1.,p[3],p[4]])
    length=np.array([1.,p[5],p[6]])
    m,g=matrices(theta,p[0],stiffness,drag,length)
    difference=np.eye(3)-np.eye(3,k=-1)
    elastic=difference.T@np.diag(stiffness)@difference
    rest=np.array([0.,p[7],p[7]+p[8]])
    return np.linalg.solve(m,g+elastic@rest)


def write(path,record):
    path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    kernel=Path(__file__).with_name('follower_curved.cl');kernel_text=kernel.read_text()
    device=cl.get_platforms()[0].get_devices()[0]
    context=cl.Context([device]);queue=cl.CommandQueue(context)
    program=cl.Program(context,kernel_text).build(options=[])
    advance=cl.Kernel(program,'advance');evaluate=cl.Kernel(program,'evaluate')
    u=qmc.Sobol(9,scramble=True,seed=17346).random_base2(args.power)
    lo=np.array([.6,.1,.1,.2,.2,.5,.5]);hi=np.array([60.,10.,10.,5.,5.,2.,2.])
    params=np.column_stack((np.exp(np.log(lo)+u[:,:7]*np.log(hi/lo)),2.8*u[:,7:]-1.4))
    initial=np.zeros((len(u),6));initial[:,1]=params[:,7];initial[:,2]=params[:,7]+params[:,8]
    initial[:,:3]+=[.1,.02,-.03];initial[:,3]=1.
    states=initial.copy();sums=np.zeros(len(u));status=np.zeros(len(u),np.int32)
    mf=cl.mem_flags
    sb=cl.Buffer(context,mf.READ_WRITE|mf.COPY_HOST_PTR,hostbuf=states)
    pb=cl.Buffer(context,mf.READ_ONLY|mf.COPY_HOST_PTR,hostbuf=params)
    lb=cl.Buffer(context,mf.READ_WRITE|mf.COPY_HOST_PTR,hostbuf=sums)
    fb=cl.Buffer(context,mf.READ_WRITE|mf.COPY_HOST_PTR,hostbuf=status)
    vb=cl.Buffer(context,mf.WRITE_ONLY,states.nbytes)
    values=np.empty_like(states)
    evaluate(queue,(len(u),),None,sb,pb,vb);cl.enqueue_copy(queue,values,vb).wait()
    worst=0.
    for i in range(len(u)):
        f=field(states[i,:3],params[i])
        jv=field(states[i,:3]+1e-25j*states[i,3:],params[i]).imag/1e-25
        worst=max(worst,float(np.max(abs(values[i]-np.r_[f,jv]))/(1+np.max(abs(np.r_[f,jv])))))
    assert worst<1e-11
    np.savez_compressed(out/'declared-inputs.npz',parameters=params,initial=initial)
    protocol=dict(seed=17346, count=len(u), parameters=['force','k2','k3','drag2','drag3','length2','length3','rest_turn2','rest_turn3'],
        logarithmic_bounds=[lo.tolist(),hi.tolist()],rest_turn_bounds=[-1.4,1.4],step=args.step,warm=args.warm,duration=args.duration,
        sampler='Scrambled Sobol, scipy',method='FP64 RK4 state and analytic directional variation; normalization every .1 tau',
        device=device.name,driver=device.driver_version,field_relative_error=worst,
        kernel_sha256=hashlib.sha256(kernel.read_bytes()).hexdigest(),
        model='Three planar rigid rods, linear torsion springs with rest curvature, local isotropic bead drag, tangential compressive tip force.',
        scope='Declared extension of De Canio2017 two-link model, not equations claimed verbatim from that paper.',
        limitations='No contact or nonlocal hydrodynamics. Fixed-step pilot is not evidence of novelty or attraction.')
    write(out/'protocol.json',protocol)
    for name in ['follower_curved.cl','follower_link.py','search_curved_follower.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    assert abs(args.warm/25-round(args.warm/25))<1e-10
    assert abs(args.duration/100-round(args.duration/100))<1e-10
    steps=int(round(25/args.step));qr=int(round(.1/args.step));assert steps%qr==0
    started=time.perf_counter();elapsed=0.;previous=np.zeros(len(u));windows=[]
    total=args.warm+args.duration
    while elapsed<total-1e-10:
        advance(queue,(len(u),),None,sb,pb,lb,fb,np.float64(args.step),np.int32(steps),np.int32(qr)).wait()
        elapsed+=25.
        if elapsed==args.warm:
            sums[:]=0.;cl.enqueue_copy(queue,lb,sums).wait()
        if elapsed>args.warm and abs((elapsed-args.warm)%(args.duration/4))<1e-9:
            cl.enqueue_copy(queue,states,sb);cl.enqueue_copy(queue,sums,lb);cl.enqueue_copy(queue,status,fb).wait()
            windows.append((sums-previous)/(args.duration/4));previous=sums.copy()
            np.savez_compressed(out/f'checkpoint-{int(elapsed):06}.npz',states=states,sums=sums,status=status,windows=windows,elapsed=elapsed)
        if elapsed%100==0:
            print('tau',elapsed,'wall_seconds',round(time.perf_counter()-started,1),flush=True)
    rates=sums/args.duration;windows=np.array(windows).T
    keep=(status==0)&np.isfinite(rates)&(windows.min(axis=1)>.01)
    np.savez_compressed(out/'complete.npz',parameters=params,initial=initial,states=states,rate=rates,windows=windows,status=status,retained=keep)
    write(out/'analysis.json',dict(completed=int((status==0).sum()),failed=int((status!=0).sum()),
        retained_ids=np.flatnonzero(keep).tolist(),retention='All four observed leading-rate windows > .01 per tau and no numerical failure.',
        retained_rates=rates[keep].tolist(),wall_seconds=time.perf_counter()-started,
        warning='No physical contact test or independent trajectory refinement performed yet.'))
    print('retained',int(keep.sum()),'of',len(u),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--power',type=int,default=9);p.add_argument('--step',type=float,default=.001)
    p.add_argument('--warm',type=float,default=100.);p.add_argument('--duration',type=float,default=400.)
    main(p.parse_args())

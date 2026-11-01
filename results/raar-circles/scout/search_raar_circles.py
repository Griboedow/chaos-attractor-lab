"""Existing Luke RAAR map on two spatial circle constraints.

This is a geometric feasibility benchmark, not a new algorithm, a model of
physical moving rings, or a priority claim. State is the actual R^3 iterate.
Parameters: beta, second radius, plane tilt, second center x,y,z.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.stats import qmc


@njit(cache=True)
def project(x,center,normal,radius):
    delta=x-center;v=delta-normal*np.dot(delta,normal);length=np.sqrt(np.dot(v,v))
    if length<1e-12:return center.copy(),np.zeros((3,3)),length,1
    u=v/length
    return center+radius*u,(radius/length)*(np.eye(3)-np.outer(normal,normal)-np.outer(u,u)),length,0


@njit(cache=True)
def step(x,params):
    beta,radius,tilt=params[:3];center=params[3:];normal=np.array([0.,-np.sin(tilt),np.cos(tilt)])
    b,jb,lb,status=project(x,center,normal,radius)
    if status:return x,np.eye(3),lb,1
    a,ja,la,status=project(2*b-x,np.zeros(3),np.array([0.,0.,1.]),1.)
    if status:return x,np.eye(3),min(la,lb),1
    result=beta*x+beta*a+(1-2*beta)*b
    jac=beta*np.eye(3)+beta*ja@(2*jb-np.eye(3))+(1-2*beta)*jb
    return result,jac,min(la,lb),0


@njit(cache=True,parallel=True)
def screen(params,warm,steps):
    n=len(params);rates=np.full((n,10),np.nan);bounds=np.empty((n,2,3));ends=np.empty((n,3));tails=np.full((n,2048,3),np.nan)
    counts=np.zeros((n,2),dtype=np.int64);minaxis=np.full(n,np.inf);maxbound=np.zeros(n)
    for k in prange(n):
        p=params[k];s=np.array([.317,-.211,.419]);v=np.array([1.,.371,-.217]);v/=np.linalg.norm(v)
        beta=p[0];constant=beta+abs(1-2*beta)*(np.linalg.norm(p[3:])+p[1]);bound=max(np.linalg.norm(s),constant/(1-beta))
        rr=np.zeros(10);lo=np.full(3,np.inf);hi=-lo
        for t in range(warm+steps):
            s,j,dist,status=step(s,p);minaxis[k]=min(minaxis[k],dist)
            if status:counts[k,0]=1;break
            if not np.isfinite(s).all():counts[k,0]=2;break
            maxbound[k]=max(maxbound[k],np.linalg.norm(s)/bound)
            v=j@v;length=np.linalg.norm(v)
            if not np.isfinite(length) or length<1e-290:counts[k,0]=3;break
            v/=length
            if t>=warm:
                rr[(t-warm)//(steps//10)]+=np.log(length);lo=np.minimum(lo,s);hi=np.maximum(hi,s);counts[k,1]+=1
            if t>=warm+steps-2048:tails[k,t-(warm+steps-2048)]=s
        if counts[k,0]==0:rates[k]=rr/(steps//10)
        ends[k]=s;bounds[k,0]=lo;bounds[k,1]=hi
    return rates,bounds,ends,tails,counts,minaxis,maxbound


def verify():
    rng=np.random.default_rng(2026110190);errors=[];forms=[]
    for _ in range(256):
        p=np.r_[rng.uniform(.2,.999),rng.uniform(.2,2),rng.uniform(.05,1.55),rng.uniform(-2,2,3)];x=rng.uniform(-2,2,3)
        new,j,dist,status=step(x,p);assert status==0
        fd=np.column_stack([(step(x+np.eye(3)[i]*1e-6,p)[0]-step(x-np.eye(3)[i]*1e-6,p)[0])/2e-6 for i in range(3)])
        errors.append(float(np.max(abs(j-fd))/max(1,np.max(abs(j)))))
        n=np.array([0.,-np.sin(p[2]),np.cos(p[2])]);b=project(x,p[3:],n,p[1])[0];rb=2*b-x
        a=project(rb,np.zeros(3),np.array([0.,0.,1.]),1.)[0];ra=2*a-rb
        forms.append(float(np.max(abs(new-(p[0]*(ra+x)/2+(1-p[0])*b)))))
    assert max(errors)<1e-6 and max(forms)<1e-13
    return dict(checks=256,max_relative_jacobian_difference=max(errors),max_reflector_formula_residual=max(forms))


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(10);checks=verify()
    u=qmc.Sobol(6,scramble=True,seed=2026110191).random_base2(args.power)
    params=np.column_stack([.2+.799*u[:,0],.2+1.8*u[:,1],.05+(np.pi/2-.05)*u[:,2],4*u[:,3:]-2])
    params=np.vstack([params,[[.8,1,0,0,0,0],[.8,1,0,3,0,0],[.8,1,np.pi/2,0,0,0]]])
    protocol=dict(count=len(params),seed=2026110191,warm=10000,observed=30000,windows=10,initial=[.317,-.211,.419],parameters=['beta','radius_B','tilt_B','center_B_x','center_B_y','center_B_z'],checks=checks,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),map='beta/2*(R_A R_B+I)+(1-beta)P_B; circles, not closed disks',primary='https://arxiv.org/pdf/math/0405208, equations (10)-(11),(14)',scope='Known RAAR algorithm, explicit geometric benchmark. Projection-axis distance below 1e-12 stops the numerical run; no arbitrary smooth extension. No novelty or asymptotic-chaos claim.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    begin=time.time();parts=[[] for _ in range(7)]
    for offset in range(0,len(params),128):
        d=screen(params[offset:offset+128],10000,30000)
        for a,b in zip(parts,d):a.append(b)
        np.savez_compressed(args.out/f'block-{offset:04}.npz',params=params[offset:offset+128],rates=d[0],bounds=d[1],ends=d[2],tails=d[3],counters=d[4],minaxis=d[5],bound_ratio=d[6])
        print(json.dumps(dict(completed=min(offset+128,len(params)),total=len(params),seconds=time.time()-begin)),flush=True)
    rates,bounds,ends,tails,counters,minaxis,bound_ratio=[np.concatenate(a) for a in parts]
    flags=(counters[:,0]==0)&(rates.min(1)>.001);assert bound_ratio.max()<1+1e-12
    np.savez_compressed(args.out/'screen.npz',params=params,rates=rates,bounds=bounds,ends=ends,tails=tails,counters=counters,minaxis=minaxis,bound_ratio=bound_ratio,flags=flags)
    report=dict(protocol,seconds=time.time()-begin,complete=int((counters[:,0]==0).sum()),flags=int(flags.sum()),max_bound_ratio=float(bound_ratio.max()))
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--power',type=int,default=11);main(p.parse_args())

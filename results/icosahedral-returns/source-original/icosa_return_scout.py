"""Finite Poincare-return scout near the numerical triangle circuit."""
from pathlib import Path
import argparse,json,time
import numpy as np
from numba import njit,prange,set_num_threads

ap=argparse.ArgumentParser();ap.add_argument('--name',default='icosa-return-scout');ap.add_argument('--count',type=int,default=513);ap.add_argument('--width',type=float,default=.015);ap.add_argument('--step',type=float,default=.005);ap.add_argument('--horizon',type=float,default=100.);args=ap.parse_args()
base=Path(__file__).parent;out=base/args.name;out.mkdir(exist_ok=False);set_num_threads(3)
d=np.load(base/'icosa-v1/input.npz');c=d['centers'][177];u=d['axes'][177];s=d['sigma'][177];gamma=d['gamma'][177]
row=json.loads((base/'icosa-cycle-refinement/report.json').read_text())['rows'][0];x0=np.array(row['point']);M=np.load(base/'icosa-cycle-refinement/cycle-0.npz')['full_monodromy']
@njit
def f(x):
    v=-gamma*x
    for j in range(60):
        r=x-c[j];a=u[j];w=np.exp(-.5*np.dot(r,r)/(s*s))
        v[0]+=w*(a[1]*r[2]-a[2]*r[1]);v[1]+=w*(a[2]*r[0]-a[0]*r[2]);v[2]+=w*(a[0]*r[1]-a[1]*r[0])
    return v
@njit
def step(x,h):
    a=f(x);b=f(x+.5*h*a);e=f(x+.5*h*b);z=f(x+h*e);return x+h/6*(a+2*b+2*e+z)
n=f(x0);n/=np.linalg.norm(n);v=np.eye(3)[np.argmin(abs(n))];v-=n*(v@n);v/=np.linalg.norm(v);Q=np.column_stack((v,np.cross(n,v)))
P=(np.eye(3)-np.outer(f(x0),n)/(n@f(x0)))@M;ev,V=np.linalg.eig(Q.T@P@Q);order=np.argsort(abs(ev))[::-1];B=Q@V[:,order].real;B/=np.linalg.norm(B,axis=0);inv=np.linalg.pinv(B)
assert np.max(abs(n@B))<1e-12 and abs(ev[order[0]]+456.17226357)<1e-5
parameters=np.linspace(-args.width,args.width,args.count);starts=x0+parameters[:,None]*B[:,0];maxevents=30
@njit(parallel=True)
def scan(starts,h,horizon):
    events=np.full((len(starts),maxevents,4),np.nan);counts=np.zeros(len(starts),np.int64);ends=np.zeros_like(starts)
    for i in prange(len(starts)):
        x=starts[i].copy();old=(x-x0)@n;t=0.;record=0
        for k in range(int(np.ceil(horizon/h))):
            y=step(x,h);new=(y-x0)@n
            if old<0 and new>=0 and t>1.:
                # Refine the actual RK4 segment rather than interpolate a chord.
                low=0.;high=h
                for _ in range(28):
                    mid=.5*(low+high);z=step(x,mid)
                    if (z-x0)@n<0:low=mid
                    else:high=mid
                z=step(x,.5*(low+high));uv=inv@(z-x0)
                if np.max(np.abs(uv))<.15 and n@f(z)>.05:
                    if record<maxevents:events[i,record,0]=t+.5*(low+high);events[i,record,1:3]=uv;events[i,record,3]=n@f(z)
                    record+=1
                    if record==maxevents:break
            x=y;old=new;t+=h
        counts[i]=record;ends[i]=x
    return events,counts,ends
begin=time.time();events,counts,ends=scan(starts,args.step,args.horizon)
np.savez_compressed(out/'data.npz',parameters=parameters,starts=starts,events=events,event_counts=counts,ends=ends,point=x0,normal=n,basis=B,inverse_basis=inv,section_derivative=inv@P@B)
report={'arguments':vars(args),'seconds':time.time()-begin,'trajectories':len(starts),'with_local_return':int((counts>0).sum()),'total_recorded_crossings':int(counts.sum()),'numerical_local_multipliers':ev[order].real.tolist(),'section_normal':n.tolist(),'section_basis':B.tolist(),'section_patch_half_width':.15,'minimum_crossing_speed':.05,'limits':'Finite RK4 scout near an approximate periodic orbit. Sampled returns and derivatives do not prove a horseshoe, invariant domain, periodic existence or chaos.'}
(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)

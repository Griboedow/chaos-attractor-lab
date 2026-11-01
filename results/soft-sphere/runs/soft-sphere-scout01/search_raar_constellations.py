"""Known RAAR iteration: unit circle A and a finite spatial constellation B.

The state is the unmodified R3 iterate. B contains polyhedron vertices, not
their convex hull. Nearest-point ties stop a scout; no smoothing is introduced.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from scipy.spatial.transform import Rotation
from numba import njit,prange,set_num_threads


def vertices(kind):
    if kind==0:a=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1.]])
    elif kind==1:a=np.vstack([np.eye(3),-np.eye(3)])
    elif kind==2:a=np.array(list(itertools.product([-1.,1.],repeat=3)))
    else:
        phi=(1+np.sqrt(5))/2
        a=np.array([[0,s,t*phi] for s in [-1,1] for t in [-1,1]]+[[t*phi,0,s] for s in [-1,1] for t in [-1,1]]+[[s,t*phi,0] for s in [-1,1] for t in [-1,1]])
    return a/np.linalg.norm(a,axis=1)[:,None]


@njit(cache=True)
def step(x,beta,points):
    closest=-1;low=np.inf;second=np.inf
    for k in range(len(points)):
        d=np.sum((x-points[k])**2)
        if d<low:second=low;low=d;closest=k
        elif d<second:second=d
    gap=second-low
    if gap<1e-12*max(1.,second):return x,np.eye(3),gap,0.,1
    b=points[closest];v=2*b-x;r=np.hypot(v[0],v[1])
    if r<1e-12:return x,np.eye(3),gap,r,2
    u=np.array([v[0]/r,v[1]/r,0.]);t=np.array([-u[1],u[0],0.])
    return beta*x+beta*u+(1-2*beta)*b,beta*(np.eye(3)-np.outer(t,t)/r),gap,r,0


@njit(cache=True,parallel=True)
def screen(params,points,counts,warm,observed):
    n=len(params);rates=np.full((n,10),np.nan);tails=np.full((n,2048,3),np.nan);ends=np.zeros((n,3))
    status=np.zeros((n,2),np.int64);minima=np.full((n,2),np.inf);ratios=np.zeros(n)
    for k in prange(n):
        beta=params[k,0];p=points[k,:counts[k]];x=np.array([.317,-.211,.419]);w=np.array([1.,.371,-.217]);w/=np.linalg.norm(w)
        radius=0.
        for v in p:radius=max(radius,np.linalg.norm(v))
        bound=max(np.linalg.norm(x),(beta+abs(1-2*beta)*radius)/(1-beta));sums=np.zeros(10)
        for it in range(warm+observed):
            x,j,gap,axis,flag=step(x,beta,p);minima[k,0]=min(minima[k,0],gap);minima[k,1]=min(minima[k,1],axis)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            w=j@w;length=np.linalg.norm(w)
            if length<1e-280 or not np.isfinite(length):status[k,0]=4;break
            w/=length;ratios[k]=max(ratios[k],np.linalg.norm(x)/bound)
            if it>=warm:sums[(it-warm)//(observed//10)]+=np.log(length);status[k,1]+=1
            if it>=warm+observed-2048:tails[k,it-(warm+observed-2048)]=x
        if status[k,0]==0:rates[k]=sums/(observed//10)
        ends[k]=x
    return rates,tails,ends,status,minima,ratios


def independent(x,beta,points):
    b=points[np.argmin(np.linalg.norm(points-x,axis=1))];rb=2*b-x;a=np.r_[rb[:2]/np.linalg.norm(rb[:2]),0.]
    return beta/2*(2*a-rb+x)+(1-beta)*b


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    u=qmc.Sobol(8,scramble=True,seed=2026110115).random_base2(args.power);params=np.column_stack([.6+.399*u[:,0],.2+2*u[:,1],1.6*u[:,2:5]-.8,u[:,5:]])
    # Every geometric draw is evaluated with each of four standard point sets.
    params=np.repeat(params,4,axis=0);kinds=np.tile(np.arange(4),len(u));points=np.zeros((len(params),12,3));counts=np.zeros(len(params),int)
    for k,(p,kind) in enumerate(zip(params,kinds)):
        a,b,c=p[5:];quat=np.array([np.sqrt(1-a)*np.sin(2*np.pi*b),np.sqrt(1-a)*np.cos(2*np.pi*b),np.sqrt(a)*np.sin(2*np.pi*c),np.sqrt(a)*np.cos(2*np.pi*c)])
        v=vertices(kind);counts[k]=len(v);points[k,:len(v)]=p[1]*Rotation.from_quat(quat).apply(v)+p[2:5]
    rng=np.random.default_rng(2026110116);forms=[];diffs=[]
    for k in range(64):
        x=rng.uniform(-2,2,3);p=points[k,:counts[k]];beta=params[k,0];y,j,_,_,flag=step(x,beta,p);assert not flag
        forms.append(float(np.max(abs(y-independent(x,beta,p)))))
        fd=np.column_stack([(independent(x+np.eye(3)[i]*1e-6,beta,p)-independent(x-np.eye(3)[i]*1e-6,beta,p))/2e-6 for i in range(3)])
        diffs.append(float(np.max(abs(j-fd))/max(1.,np.max(abs(j)))))
    assert max(forms)<1e-12 and max(diffs)<1e-5
    protocol=dict(seed=2026110115,count=len(params),kinds=['tetrahedron vertices','octahedron vertices','cube vertices','icosahedron vertices'],warm=5000,observed=15000,windows=10,initial=[.317,-.211,.419],checks=dict(n=64,nested_formula=max(forms),relative_fd_jacobian=max(diffs)),map='beta/2*(R_A R_B+I)+(1-beta)P_B',source='https://arxiv.org/abs/math/0405208',scope='Known algorithm. Piecewise smooth map with Voronoi jumps. Local derivative growth does not certify chaos. No novelty claim.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());begin=time.time();parts=[[] for _ in range(6)]
    for offset in range(0,len(params),128):
        d=screen(params[offset:offset+128],points[offset:offset+128],counts[offset:offset+128],5000,15000)
        for a,b in zip(parts,d):a.append(b)
        np.savez_compressed(args.out/f'block-{offset:04}.npz',params=params[offset:offset+128],points=points[offset:offset+128],counts=counts[offset:offset+128],kinds=kinds[offset:offset+128],rates=d[0],tails=d[1],ends=d[2],status=d[3],minima=d[4],boundratios=d[5])
        print(json.dumps(dict(completed=min(offset+128,len(params)),total=len(params),seconds=time.time()-begin)),flush=True)
    rates,tails,ends,status,minima,ratios=[np.concatenate(v) for v in parts];flags=(status[:,0]==0)&(rates.min(1)>.005);assert ratios.max()<1+1e-12
    np.savez_compressed(args.out/'screen.npz',params=params,points=points,counts=counts,kinds=kinds,rates=rates,tails=tails,ends=ends,status=status,minima=minima,boundratios=ratios,flags=flags)
    report=dict(protocol,complete=int((status[:,0]==0).sum()),flagged=int(flags.sum()),flagged_ids=np.flatnonzero(flags).tolist(),max_boundratio=float(ratios.max()),seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--power',type=int,default=9);main(p.parse_args())

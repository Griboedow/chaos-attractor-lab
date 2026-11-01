"""Paired RAAR survey: replace the circle by the unit sphere in native R3.

All parameter and finite-set bytes come from the declared circle survey.
Voronoi ties and the sphere projection singularity stop the calculation.
"""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from numba import njit, prange, set_num_threads


@njit(cache=True)
def step(x, beta, points):
    closest=-1; low=np.inf; second=np.inf
    for k in range(len(points)):
        d=np.sum((x-points[k])**2)
        if d<low: second=low; low=d; closest=k
        elif d<second: second=d
    gap=second-low
    if gap<1e-12*max(1.,second): return x,np.eye(3),gap,0.,1
    b=points[closest]; v=2*b-x; r=np.linalg.norm(v)
    if r<1e-12: return x,np.eye(3),gap,r,2
    u=v/r
    return beta*x+beta*u+(1-2*beta)*b,beta*(np.eye(3)-(np.eye(3)-np.outer(u,u))/r),gap,r,0


@njit(cache=True,parallel=True)
def screen(params,points,counts,warm,observed):
    n=len(params); rates=np.full((n,10),np.nan); tails=np.full((n,2048,3),np.nan)
    status=np.zeros((n,2),np.int64); minima=np.full((n,2),np.inf); ratios=np.zeros(n)
    for k in prange(n):
        beta=params[k,0]; p=points[k,:counts[k]]; x=np.array([.317,-.211,.419])
        w=np.array([1.,.371,-.217]); w/=np.linalg.norm(w); radius=0.
        for v in p: radius=max(radius,np.linalg.norm(v))
        bound=max(np.linalg.norm(x),(beta+abs(1-2*beta)*radius)/(1-beta)); sums=np.zeros(10)
        for it in range(warm+observed):
            x,j,gap,r,flag=step(x,beta,p); minima[k,0]=min(minima[k,0],gap); minima[k,1]=min(minima[k,1],r)
            if flag or not np.isfinite(x).all(): status[k,0]=flag if flag else 3; break
            w=j@w; length=np.linalg.norm(w)
            if length<1e-280 or not np.isfinite(length): status[k,0]=4; break
            w/=length; ratios[k]=max(ratios[k],np.linalg.norm(x)/bound)
            if it>=warm: sums[(it-warm)//(observed//10)]+=np.log(length); status[k,1]+=1
            if it>=warm+observed-2048: tails[k,it-(warm+observed-2048)]=x
        if status[k,0]==0: rates[k]=sums/(observed//10)
    return rates,tails,status,minima,ratios


def independent(x,beta,points):
    b=points[np.argmin(np.linalg.norm(points-x,axis=1))]; rb=2*b-x; a=rb/np.linalg.norm(rb)
    return beta/2*(2*a-rb+x)+(1-beta)*b


def main(args):
    args.out.mkdir(parents=True,exist_ok=False); set_num_threads(4)
    with np.load(args.source) as f:
        params=f['params']; points=f['points']; counts=f['counts']; kinds=f['kinds']; circle_flags=f['flags']
    rng=np.random.default_rng(2026110118); forms=[]; diffs=[]
    for k in range(64):
        x=rng.uniform(-2,2,3); p=points[k,:counts[k]]; beta=params[k,0]; y,j,_,_,flag=step(x,beta,p); assert not flag
        forms.append(float(np.max(abs(y-independent(x,beta,p)))))
        fd=np.column_stack([(independent(x+np.eye(3)[i]*1e-6,beta,p)-independent(x-np.eye(3)[i]*1e-6,beta,p))/2e-6 for i in range(3)])
        diffs.append(float(np.max(abs(j-fd))/max(1.,np.max(abs(j)))))
    assert max(forms)<1e-12 and max(diffs)<1e-5
    protocol=dict(input=str(args.source),input_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),count=len(params),warm=5000,observed=15000,windows=10,initial=[.317,-.211,.419],checks=dict(n=64,nested_formula=max(forms),relative_fd_jacobian=max(diffs)),map='beta/2*(R_A R_B+I)+(1-beta)P_B; A=unit sphere; B=finite vertices',scope='Known RAAR formula with a changed constraint. Paired input parameters are identical to the circle survey. Discontinuous map; local derivative growth does not certify chaos or novelty.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes()); begin=time.time(); parts=[[] for _ in range(5)]
    for offset in range(0,len(params),128):
        s=slice(offset,offset+128); d=screen(params[s],points[s],counts[s],5000,15000)
        for a,b in zip(parts,d): a.append(b)
        np.savez_compressed(args.out/f'block-{offset:04}.npz',params=params[s],points=points[s],counts=counts[s],kinds=kinds[s],rates=d[0],tails=d[1],status=d[2],minima=d[3],boundratios=d[4])
        print(json.dumps(dict(completed=min(offset+128,len(params)),total=len(params),seconds=time.time()-begin)),flush=True)
    rates,tails,status,minima,ratios=[np.concatenate(v) for v in parts]
    flags=(status[:,0]==0)&(rates.min(1)>.005); assert ratios.max()<1+1e-12
    np.savez_compressed(args.out/'screen.npz',params=params,points=points,counts=counts,kinds=kinds,rates=rates,tails=tails,status=status,minima=minima,boundratios=ratios,flags=flags,circle_flags=circle_flags)
    report=dict(protocol,complete=int((status[:,0]==0).sum()),flagged=int(flags.sum()),flagged_ids=np.flatnonzero(flags).tolist(),both_flagged_ids=np.flatnonzero(flags&circle_flags).tolist(),max_boundratio=float(ratios.max()),seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n'); print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--source',type=Path,required=True); p.add_argument('--out',type=Path,required=True); main(p.parse_args())

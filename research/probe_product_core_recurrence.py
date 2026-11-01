"""Streaming recurrence search towards the certified core, in native coordinates."""
import argparse,hashlib,json,time
from pathlib import Path
from fractions import Fraction
import numpy as np
from numba import njit,prange,set_num_threads


@njit(cache=True)
def value(x,beta,tau,levels):
    m=np.empty(3)
    for j in range(3):
        lo,hi=levels[j,:2];a=(hi-lo)/2;c=(hi+lo)/2
        m[j]=c+a*np.tanh(2*a*(x[j]-c)/tau)
    v=2*m-x;r=np.linalg.norm(v)
    assert r>1e-12
    return beta*x+beta*v/r+(1-2*beta)*m,r


@njit(cache=True)
def run(start,beta,tau,levels,targets,inverse0,warm,count,block,radii):
    x=start.copy();best=np.full(len(targets),np.inf);best_states=np.zeros((len(targets),3))
    best_it=np.zeros(len(targets),np.int64);counts=np.zeros((count//block,len(radii)),np.int64)
    history=np.empty((count//block,len(targets)));minimum=np.inf;window_hits=0
    for it in range(warm+count):
        x,r=value(x,beta,tau,levels);minimum=min(minimum,r)
        if it<warm:continue
        n=it-warm;ib=n//block
        for k in range(len(targets)):
            d=np.sum((x-targets[k])**2)
            if d<best[k]:best[k]=d;best_states[k]=x;best_it[k]=n+1
            if k==0:
                for j in range(len(radii)):
                    if d<radii[j]*radii[j]:counts[ib,j]+=1
                if d<1e-8:
                    local=inverse0@(x-targets[0])
                    if np.max(np.abs(local))<=1:window_hits+=1
        if (n+1)%block==0:history[ib]=np.sqrt(best)
    return x,history,best_states,best_it,counts,minimum,window_hits


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    raw=json.loads(args.input.read_text());cycles=np.load(args.cycle)
    beta,tau=cycles['params'];levels=cycles['levels']
    # Four phases of cycle 14 and the other 22 covering-window centers.
    targets=np.concatenate([cycles['orbit'],np.array(raw['centers'])[1:]])
    inverse0=np.array([[float(Fraction(v)) for v in row] for row in raw['exact_windows'][0]['inverse']])
    radii=np.array([.005,.002,.001,.0005,.0001,.00001]);starts=np.array([[.317,-.211,.419],[-.42,.18,-.31]])
    reports=[];payload={}
    for k,start in enumerate(starts):
        end,history,best,when,counts,minimum,hits=run(start,beta,tau,levels,targets,inverse0,20000,args.count,args.block,radii)
        payload.update({f'history_{k}':history,f'closest_states_{k}':best,f'closest_iterations_{k}':when,f'counts_{k}':counts,f'last_{k}':end})
        reports.append(dict(start=start.tolist(),closest_distances=history[-1].tolist(),closest_iterations=when.tolist(),
                            counts_by_radius=counts.sum(0).tolist(),minimum_normalization=float(minimum),window0_hits=int(hits)))
        print(json.dumps(dict(completed=k+1,seconds=time.time()-begin,report=reports[-1])),flush=True)
    np.savez_compressed(args.out/'recurrence.npz',targets=targets,starts=starts,radii=radii,**payload)
    report=dict(status='completed',sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.input,args.cycle]},
        warm=20000,count_per_start=args.count,block=args.block,radii=radii.tolist(),target_order='Cycle-14 phases 0..3, then windows 1..22.',
        formula='Binary tanh posterior with float64 half-sums/differences; independent floating realization, not an exact-orbit certificate.',
        starts=reports,seconds=time.time()-begin,
        scope='Finite recurrence search. A plateau in nearest distance is evidence against identification with the observed recurrent support, not proof of a positive invariant gap.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','cycle','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--count',type=int,default=50000000);p.add_argument('--block',type=int,default=1000000)
    main(p.parse_args())

"""Declared starts and parameter neighbours for six selected spatial-circle RAAR cases."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.stats import qmc
from search_raar_circles import step
from follow_raar_circles import orthogonalize

@njit(cache=True)
def one(params,start,warm,steps):
    state=start.copy();Q,_,_=orthogonalize(np.array([[1.,.2,-.3],[.4,1.,.15],[-.25,.35,1.]]))
    rates=np.zeros((10,3));bounds=np.array([[np.inf,np.inf,np.inf],[-np.inf,-np.inf,-np.inf]])
    tail=np.full((2048,3),np.nan);count=np.zeros(3,dtype=np.int64)
    for n in range(warm+steps):
        if not np.all(np.isfinite(state)) or np.max(np.abs(state))>1e10:count[0]=2;break
        state,J,axis_distance,status=step(state,params);reject=0
        if status:count[0]=1;break
        Q,r,status=orthogonalize(J@Q)
        if status:count[0]=3;break
        if n>=warm:
            rates[(n-warm)//(steps//10)]+=r;count[1]+=1;count[2]+=reject
            bounds[0]=np.minimum(bounds[0],state);bounds[1]=np.maximum(bounds[1],state)
        if n>=warm+steps-2048:tail[n-(warm+steps-2048)]=state
    return rates/(steps//10),bounds,tail,count

@njit(cache=True,parallel=True)
def batch(params,starts,warm,steps):
    n=len(params);rates=np.empty((n,10,3));bounds=np.empty((n,2,3));tails=np.empty((n,2048,3));counts=np.empty((n,3),dtype=np.int64)
    for i in prange(n):rates[i],bounds[i],tails[i],counts[i]=one(params[i],starts[i],warm,steps)
    return rates,bounds,tails,counts

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6)
    with np.load(args.source) as f:source=f['params']
    params=[];starts=[];ids=[];groups=[];details=[]
    for ident in [75,820,1051,1311,1483,1679]:
        base=source[ident]
        def add(p,x,g,d):params.append(p);starts.append(x);ids.append(ident);groups.append(g);details.append(d)
        raw=qmc.Sobol(3,scramble=True,seed=2026110300+ident).random_base2(5)
        for z in raw:add(base,8*z-4,'start',dict(box=[-4,4]))
        neighbors=qmc.Sobol(6,scramble=True,seed=2026110400+ident).random_base2(5)
        for z in neighbors:
            p=base+.01*(2*z-1)
            assert 0<p[0]<1 and p[1]>0
            add(p,[.317,-.211,.419],'neighbor',dict(absolute_halfwidth=.01))
    params=np.array(params);starts=np.array(starts);ids=np.array(ids);groups=np.array(groups)
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),warm=args.warm,observed=args.steps,windows=10,cases=[75,820,1051,1311,1483,1679],count=len(ids),details=details,scope='Finite basin and parameter-neighbourhood evidence, not basin measures or structural-stability proof.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    for name in ['cohort_raar_circles.py','follow_raar_circles.py','search_raar_circles.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    np.savez_compressed(args.out/'inputs.npz',params=params,starts=starts,ids=ids,groups=groups)
    begin=time.time();result=batch(params,starts,args.warm,args.steps);rates,bounds,tails,counters=result
    period=np.zeros(len(ids),dtype=int)
    for i in range(len(ids)):
        if counters[i,0]:continue
        for lag in range(1,513):
            if np.max(np.abs(tails[i,lag:]-tails[i,:-lag]))<1e-9:period[i]=lag;break
    flags=(counters[:,0]==0)&(rates[:,:,0].min(1)>.001)&(period==0)
    np.savez_compressed(args.out/'cohort.npz',params=params,starts=starts,ids=ids,groups=groups,rates=rates,bounds=bounds,tails=tails,counters=counters,period=period,flags=flags)
    records=[]
    for ident in [75,820,1051,1311,1483,1679]:
        for group in np.unique(groups):
            mask=(ids==ident)&(groups==group)
            records.append(dict(id=ident,group=str(group),count=int(mask.sum()),completed=int(np.sum(mask&(counters[:,0]==0))),positive=int(np.sum(mask&flags)),periodic=int(np.sum(mask&(period>0))),two_positive=int(np.sum(mask&flags&(rates[:,:,1].min(1)>0)))))
    report=dict(protocol,seconds=time.time()-begin,groups=records);(args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--warm',type=int,default=100000);p.add_argument('--steps',type=int,default=300000);main(p.parse_args())

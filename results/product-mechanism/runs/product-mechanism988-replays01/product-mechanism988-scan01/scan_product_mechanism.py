"""Declared two-parameter slices, three starts each; finite observations only."""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from numba import njit, prange, set_num_threads
from search_product_alphabet import step


@njit(cache=True, parallel=True)
def survey(params, levels, starts, warm, count):
    n=len(params);rates=np.zeros((n,len(starts),10));tails=np.zeros((n,len(starts),2048,3))
    minima=np.full((n,len(starts)),np.inf);flags=np.zeros((n,len(starts)),np.int64)
    active=np.zeros((n,len(starts),4),np.int64)
    for k in prange(n):
        beta,tau=params[k]
        for s in range(len(starts)):
            x=starts[s].copy();w=np.array([1.,.371,-.217]);w/=np.linalg.norm(w)
            for it in range(warm+count):
                y,J,r,_,_,flag=step(x,beta,tau,levels,2)
                if flag or not np.isfinite(y).all():flags[k,s]=1;break
                minima[k,s]=min(minima[k,s],r)
                w=J@w;length=np.linalg.norm(w)
                if length<1e-280 or not np.isfinite(length):flags[k,s]=2;break
                w/=length
                if it>=warm:
                    rates[k,s,(it-warm)//(count//10)]+=np.log(length)/(count//10)
                    switching=0
                    for j in range(3):
                        a=(levels[j,1]-levels[j,0])/2;c=(levels[j,1]+levels[j,0])/2
                        e=np.exp(-abs(4*a*(x[j]-c)/tau));d=8*a*a/tau*e/(1+e)**2
                        if d>1:switching+=1
                    active[k,s,switching]+=1
                x=y
                if it>=warm+count-2048:tails[k,s,it-(warm+count-2048)]=x
    return rates,tails,minima,flags,active


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time();set_num_threads(4)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    with np.load(args.source) as f:beta,tau=f['params'];levels=f['levels'][:,:2]
    starts=np.array([[.317,-.211,.419],[-.42,.18,-.31],[.11,.33,-.27]])
    beta_grid=np.linspace(.88,.985,args.size);tau_grid=np.geomspace(.004,.06,args.size)
    # The exact principal case is retained in each slice, without replacing a grid point.
    bg=np.sort(np.r_[beta_grid,beta]);tg=np.sort(np.r_[tau_grid,tau])
    params=np.vstack([np.column_stack([bg,np.full(len(bg),tau)]),np.column_stack([np.full(len(tg),beta),tg])])
    protocol=dict(source_sha256=sha(args.source),beta=float(beta),tau=float(tau),levels=levels.tolist(),starts=starts.tolist(),
                  parameters=params.tolist(),slice_lengths=[len(bg),len(tg)],warm=args.warm,observed=args.count,
                  scope='Independent starts at every parameter, not continuation. Finite leading growth and last 2048 states; no bifurcation or basin theorem.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rates,tails,minima,flags,active=survey(params,levels,starts,args.warm,args.count)
    periods=np.zeros(flags.shape,np.int64)
    for k in range(len(params)):
        for s in range(len(starts)):
            if flags[k,s]:continue
            for p in range(1,513):
                if np.max(abs(tails[k,s,-512:]-tails[k,s,-512-p:-p]))<1e-9:periods[k,s]=p;break
    np.savez_compressed(args.out/'scan.npz',params=params,levels=levels,starts=starts,rates=rates,tails=tails,minima=minima,flags=flags,active=active,periods=periods)
    rows=[]
    for k in range(len(params)):
        rows.append(dict(id=k,params=params[k].tolist(),mean_rates=rates[k].mean(1).tolist(),minimum_blocks=rates[k].min(1).tolist(),
                         flags=flags[k].tolist(),finite_periods=periods[k].tolist(),minimum_radii=minima[k].tolist(),active_counts=active[k].tolist()))
    report=dict(status='completed',protocol_sha256=sha(args.out/'protocol.json'),cases=rows,seconds=time.time()-begin,
                strictly_positive_cohorts=int((rates.min(2)>0).sum()),finite_periodic_cohorts=int((periods>0).sum()),flagged_cohorts=int((flags>0).sum()),scope=protocol['scope'])
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report[k] for k in ['status','seconds','strictly_positive_cohorts','finite_periodic_cohorts','flagged_cohorts']}))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--size',type=int,default=241);p.add_argument('--warm',type=int,default=20000);p.add_argument('--count',type=int,default=80000)
    main(p.parse_args())

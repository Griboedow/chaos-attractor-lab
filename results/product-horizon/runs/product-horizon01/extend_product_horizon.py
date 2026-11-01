"""Stream a longer numerical horizon for the two already selected product cases.

Every block retains its incoming and outgoing state/frame, finite-time spectrum,
normalization minimum and bounds. This is not an asymptotic certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit, prange, set_num_threads
from search_product_alphabet import step
from follow_spherical_constellations import qr_step


@njit(cache=True, parallel=True)
def advance(states, frames, params, levels, orders, count, keep):
    n=len(states); out=states.copy(); qs=frames.copy(); logs=np.zeros((n,3))
    minima=np.full(n,np.inf); bounds=np.full((n,2,3),np.inf); bounds[:,1,:]=-np.inf
    tails=np.empty((n,keep,3)); flags=np.zeros((n,2),np.int64)
    for k in prange(n):
        x=out[k].copy(); Q=qs[k].copy()
        for t in range(count):
            y,J,r,_,_,flag=step(x,params[k,0],params[k,1],levels[k],orders[k])
            if flag or not np.isfinite(y).all(): flags[k,0]=1;break
            Q,ell,flag=qr_step(J,Q)
            if flag: flags[k,0]=2;break
            x=y;logs[k]+=ell;minima[k]=min(minima[k],r);flags[k,1]+=1
            for j in range(3): bounds[k,0,j]=min(bounds[k,0,j],x[j]);bounds[k,1,j]=max(bounds[k,1,j],x[j])
            if t>=count-keep: tails[k,t-(count-keep)]=x
        out[k]=x;qs[k]=Q
    return out,qs,logs,minima,bounds,tails,flags


def save(p, obj): p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def main(args):
    args.out.mkdir(parents=True,exist_ok=False); start=time.time();set_num_threads(2)
    names=[Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']
    for name in names: (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    params=[];levels=[];orders=[];inputs=[]
    for ident in args.ids:
        source=args.source/f'case-{ident}-states.npz'
        with np.load(source) as f: params.append(f['params']);levels.append(f['levels']);orders.append(int(f['order']))
        inputs.append(dict(id=ident,path=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    params=np.array(params);levels=np.array(levels);orders=np.array(orders)
    x=np.tile([.317,-.211,.419],(len(inputs),1));Q=np.tile(np.eye(3),(len(inputs),1,1))
    protocol=dict(inputs=inputs,ids=args.ids,warm=args.warm,observed=args.count,block=args.block,initial=x.tolist(),params=params.tolist(),levels=levels.tolist(),orders=orders.tolist(),selection='The two already published principal portraits, chosen before this extended-horizon calculation.',scope='Fresh-start binary64 finite-time continuation. No proof of asymptotic attraction, transitivity or uniqueness.')
    save(args.out/'protocol.json',protocol)
    x,Q,_,_,_,_,flags=advance(x,Q,params,levels,orders,args.warm,0);assert not np.any(flags[:,0])
    assert args.count%args.block==0 and args.block>=100000
    rows=[]; totals=np.zeros((len(inputs),3));minima=np.full(len(inputs),np.inf)
    for block in range(args.count//args.block):
        incoming=x.copy();incomingQ=Q.copy();keep=100000 if (block+1)*args.block==args.count else 0
        x,Q,logs,r,bounds,tail,flags=advance(x,Q,params,levels,orders,args.block,keep)
        assert not np.any(flags[:,0]) and np.all(flags[:,1]==args.block),flags
        totals+=logs;minima=np.minimum(minima,r)
        row=dict(block=block,start=block*args.block,count=args.block,incoming_state=incoming.tolist(),incoming_frame=incomingQ.tolist(),outgoing_state=x.tolist(),outgoing_frame=Q.tolist(),spectrum=(logs/args.block).tolist(),minimum_radius=r.tolist(),bounds=bounds.tolist())
        rows.append(row);save(args.out/f'block-{block:03}.json',row)
        save(args.out/'progress.json',dict(completed=(block+1)*args.block,total=args.count,mean_spectrum=(totals/((block+1)*args.block)).tolist(),minimum_radius=minima.tolist(),seconds=time.time()-start))
        print(json.dumps(dict(completed=(block+1)*args.block,total=args.count,spectra=(logs/args.block).tolist(),seconds=time.time()-start)),flush=True)
    cases=[]
    for k,ident in enumerate(args.ids):
        z=tail[k];scale=max(np.linalg.norm(np.ptp(z,axis=0)),1.)
        residual=[float(np.max(np.linalg.norm(z[-2048+lag:]-z[-2048:-lag],axis=1))) for lag in range(1,513)]
        periods=[i+1 for i,v in enumerate(residual) if v<1e-8*scale]
        ev=np.linalg.eigvalsh(np.cov(z.T))
        cases.append(dict(id=ident,spectrum=(totals[k]/args.count).tolist(),minimum_leading_block=min(row['spectrum'][k][0] for row in rows),maximum_leading_block=max(row['spectrum'][k][0] for row in rows),minimum_radius=float(minima[k]),detected_period=min(periods) if periods else None,return_residuals=residual,tail_variance_fractions=(ev/ev.sum()).tolist()))
        np.savez_compressed(args.out/f'case-{ident}-tail.npz',states=z,params=params[k],levels=levels[k],order=orders[k],final_frame=Q[k])
    save(args.out/'report.json',dict(protocol,cases=cases,blocks=len(rows),seconds=time.time()-start));print(json.dumps(cases),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',default=[988,489]);p.add_argument('--count',type=int,default=100000000)
    p.add_argument('--warm',type=int,default=100000);p.add_argument('--block',type=int,default=2000000)
    main(p.parse_args())

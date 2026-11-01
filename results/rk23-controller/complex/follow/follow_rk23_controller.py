"""Fresh-start full two-dimensional spectra for every RK23 scout flag."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_rk23_controller import rk23_step

@njit(cache=True)
def one(params,warm,steps,windows,tail_count):
    s=np.array([np.log(.1),0.]);phase=0.;v=np.array([.9,.435889894]);v/=np.sqrt(np.sum(v*v))
    rates=np.zeros((windows,2));lo=np.full(2,np.inf);hi=np.full(2,-np.inf);tail=np.empty((tail_count,3));rejects=0
    for t in range(warm+steps):
        s,j,dphi,reject,status=rk23_step(s,params)
        if status or not np.all(np.isfinite(s)) or np.max(np.abs(s))>100:
            return rates,lo,hi,s,phase,tail,-1
        phase=(phase+dphi)%(2*np.pi)
        v=j@v;norm=np.sqrt(np.sum(v*v));v/=norm
        det=j[0,0]*j[1,1]-j[0,1]*j[1,0]
        if t>=warm:
            index=(t-warm)//(steps//windows);rate=np.log(norm)
            rates[index,0]+=rate;rates[index,1]+=np.log(abs(det))-rate
            lo=np.minimum(lo,s);hi=np.maximum(hi,s);rejects+=reject
        if t>=warm+steps-tail_count:tail[t-(warm+steps-tail_count)]=np.array([s[0],s[1],phase])
    rates/=steps//windows
    return rates,lo,hi,s,phase,tail,rejects

@njit(cache=True,parallel=True)
def many(params,warm,steps,windows,tail_count):
    n=len(params);rates=np.empty((n,windows,2));bounds=np.empty((n,2,2));ends=np.empty((n,3));states=np.empty((n,tail_count,3));rejects=np.empty(n,dtype=np.int64)
    for i in prange(n):
        r,lo,hi,s,phase,tail,reject=one(params[i],warm,steps,windows,tail_count)
        rates[i]=r;bounds[i,0]=lo;bounds[i,1]=hi;ends[i]=np.array([s[0],s[1],phase]);states[i]=tail;rejects[i]=reject
    return rates,bounds,ends,states,rejects

@njit(cache=True)
def cloud(params,start,count):
    s=start[:2].copy();phase=start[2];states=np.empty((count,3));logs=np.empty((count,2))
    for i in range(count):
        s,_,dphi,_,status=rk23_step(s,params)
        if status:raise ValueError('Step failed')
        phase=(phase+dphi)%(2*np.pi);amplitude=np.exp(s[1]);logs[i]=s
        states[i]=np.array([amplitude*np.cos(phase),amplitude*np.sin(phase),np.exp(s[0])])
    return states,logs

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(8)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    ids=np.flatnonzero(d['flags']);params=d['params'][ids];begin=time.time()
    rates,bounds,ends,states,rejects=many(params,100000,1000000,20,16384)
    period=np.zeros(len(ids),dtype=np.int64)
    for i,z in enumerate(states):
        if rejects[i]<0:continue
        for lag in range(1,513):
            if np.max(np.abs(z[-4096+lag:,:2]-z[-4096:-lag,:2]))<1e-9:period[i]=lag;break
    flags=(rejects>=0)&(rates[:,:,0].min(1)>.001)&(period==0)
    np.savez_compressed(args.out/'follow.npz',ids=ids,params=params,rates=rates,bounds=bounds,ends=ends,tails=states,rejects=rejects,period=period,flags=flags)
    report=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),warm=100000,observed=1000000,windows=20,
                all_scout_flags_followed=len(ids),retained=int(flags.sum()),periods={str(k):int(np.sum(period==k)) for k in np.unique(period) if k},
                method='Leading QR direction plus determinant identity for second quotient rate; independent phase symmetry adds zero.',seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    for name in ['search_rk23_controller.py','follow_rk23_controller.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

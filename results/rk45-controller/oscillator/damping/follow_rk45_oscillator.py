"""Full-spectrum fresh-start follow-up, retaining every oscillator scout flag."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_rk45_oscillator import matrices,oscillator_step

@njit(cache=True)
def orthogonalize(q):
    out=np.empty((3,3));rates=np.empty(3)
    for j in range(3):
        v=q[:,j].copy()
        for repeat in range(2):
            for k in range(j):v-=np.dot(out[:,k],v)*out[:,k]
        norm=np.sqrt(np.dot(v,v))
        if norm<1e-290:return out,rates,1
        out[:,j]=v/norm;rates[j]=np.log(norm)
    return out,rates,0

@njit(cache=True)
def one(params,warm,steps,windows):
    ms=matrices(params[0]);s=np.array([1.,0.,np.log(.1)])
    q,_,_=orthogonalize(np.array([[1.,.2,-.3],[.4,1.,.15],[-.25,.35,1.]]))
    rates=np.zeros((windows,3));detlogs=np.zeros(windows);lo=np.full(3,np.inf);hi=np.full(3,-np.inf)
    tails=np.full((2048,3),np.nan);counts=np.zeros(3,dtype=np.int64)
    for t in range(warm+steps):
        if not np.all(np.isfinite(s)) or max(abs(s[0]),abs(s[1]))>1e10 or abs(s[2])>50:counts[0]=2;break
        s,j,reject,status=oscillator_step(s,params[1],*ms)
        if status:counts[0]=1;break
        q,r,status=orthogonalize(j@q)
        if status:counts[0]=3;break
        if t>=warm:
            k=(t-warm)//(steps//windows);rates[k]+=r
            det=j[0,0]*(j[1,1]*j[2,2]-j[1,2]*j[2,1])-j[0,1]*(j[1,0]*j[2,2]-j[1,2]*j[2,0])+j[0,2]*(j[1,0]*j[2,1]-j[1,1]*j[2,0])
            detlogs[k]+=np.log(abs(det));lo=np.minimum(lo,s);hi=np.maximum(hi,s);counts[1]+=1;counts[2]+=reject
        if t>=warm+steps-2048:tails[t-(warm+steps-2048)]=s
    rates/=steps//windows;detlogs/=steps//windows
    return rates,detlogs,lo,hi,s,tails,counts

@njit(cache=True,parallel=True)
def batch(params,warm,steps,windows):
    n=len(params);rates=np.empty((n,windows,3));detlogs=np.empty((n,windows));bounds=np.empty((n,2,3));ends=np.empty((n,3));tails=np.empty((n,2048,3));counts=np.empty((n,3),dtype=np.int64)
    for i in prange(n):
        r,det,lo,hi,s,tail,count=one(params[i],warm,steps,windows)
        rates[i]=r;detlogs[i]=det;bounds[i,0]=lo;bounds[i,1]=hi;ends[i]=s;tails[i]=tail;counts[i]=count
    return rates,detlogs,bounds,ends,tails,counts

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(12)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    ids=np.flatnonzero(d['flags']) if not args.ids else np.array(args.ids,dtype=int)
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),warm=args.warm,observed=args.steps,windows=20,initialization='Fresh q=1,p=0,proposed h=.1',method='Full three-column QR with two-pass modified Gram-Schmidt',source_sha256_code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    for name in ['follow_rk45_oscillator.py','search_rk45_oscillator.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    arrays=[[] for _ in range(6)];begin=time.time()
    for offset in range(0,len(ids),24):
        which=ids[offset:offset+24];result=batch(d['params'][which],args.warm,args.steps,20)
        np.savez_compressed(args.out/f'block-{offset:04}.npz',ids=which,params=d['params'][which],rates=result[0],detlogs=result[1],bounds=result[2],ends=result[3],tails=result[4],counters=result[5])
        for a,r in zip(arrays,result):a.append(r)
        print(json.dumps(dict(completed=min(offset+24,len(ids)),total=len(ids),seconds=time.time()-begin)),flush=True)
    rates,detlogs,bounds,ends,tails,counters=[np.concatenate(a) for a in arrays];period=np.zeros(len(ids),dtype=np.int64)
    for i,z in enumerate(tails):
        if counters[i,0]:continue
        for lag in range(1,513):
            if np.max(abs(z[lag:]-z[:-lag]))<1e-9:period[i]=lag;break
    flags=(counters[:,0]==0)&(rates[:,:,0].min(1)>.001)&(period==0)
    np.savez_compressed(args.out/'follow.npz',ids=ids,params=d['params'][ids],rates=rates,detlogs=detlogs,bounds=bounds,ends=ends,tails=tails,counters=counters,period=period,flags=flags)
    report=dict(protocol,count=len(ids),retained=int(flags.sum()),two_positive=int(np.sum(flags&(rates[:,:,1].min(1)>.001))),period_counts={str(k):int(np.sum(period==k)) for k in np.unique(period) if k},seconds=time.time()-begin,
                max_window_trace_identity_error=float(np.max(abs(rates.sum(2)-detlogs))))
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+');p.add_argument('--warm',type=int,default=100000);p.add_argument('--steps',type=int,default=1000000);main(p.parse_args())

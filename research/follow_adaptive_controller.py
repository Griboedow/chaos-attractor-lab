"""Fresh-start follow-up of every PI-controller scout flag."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_adaptive_controller import advance,batch

@njit(cache=True,parallel=True)
def tails(params,policies,ends,count):
    data=np.full((len(params),count,3),np.nan)
    for i in prange(len(params)):
        z=ends[i].copy()
        for j in range(count):
            if not np.all(np.isfinite(z)) or np.max(np.abs(z))>100:
                break
            z,_,_,_=advance(z,params[i],policies[i])
            data[i,j]=z
    return data

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(12)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    ids=np.flatnonzero(d['flags']);params=d['params'][ids];policies=d['policies'][ids];initials=d['initials'][ids]
    begin=time.time();rates,bounds,ends,counters=batch(params,policies,initials,100000,500000,20)
    states=tails(params,policies,ends,4096)
    period=np.zeros(len(ids),dtype=np.int64)
    for i,z in enumerate(states):
        if counters[i,0]!=0 or not np.isfinite(z).all():continue
        for lag in range(1,513):
            if np.max(np.abs(z[lag:]-z[:-lag]))<1e-9:
                period[i]=lag;break
    flags=(counters[:,0]==0)&(rates.min(1)>.001)&(period==0)
    np.savez_compressed(args.out/'follow.npz',ids=ids,params=params,policies=policies,initials=initials,
                        rates=rates,bounds=bounds,ends=ends,counters=counters,tails=states,period=period,flags=flags)
    report={'source':str(args.source),'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),
            'source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'warm':100000,'observed':500000,'windows':20,'tail':4096,'initialization':'fresh, same as scout',
            'all_scout_flags_followed':len(ids),'retained':int(flags.sum()),'completed':int(np.sum(counters[:,0]==0)),
            'finite_periods_detected':{str(k):int(np.sum(period==k)) for k in np.unique(period) if k},
            'seconds':time.time()-begin}
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2),flush=True)
    for name in ['search_adaptive_controller.py','follow_adaptive_controller.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

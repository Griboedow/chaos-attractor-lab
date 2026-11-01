"""Declared fine slice around the observed R1051 period-doubling sequence."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from follow_raar_circles import batch


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(8)
    with np.load(args.source) as f:base=f['params'][1051]
    values=np.linspace(.9608,.9664,1025);params=np.tile(base,(len(values),1));params[:,0]=values
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),base_id=1051,beta_bounds=[.9608,.9664],count=len(params),warm=100000,observed=400000,windows=20,initial=[.317,-.211,.419],selection='Fine slice selected after the declared broad slice and periodic-orbit continuation exposed a 7/14/28/56 sequence; not an independent discovery test.',scope='Finite fresh-start outcomes. No certified bifurcation or cascade-accumulation claim.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    for name in ['scan_raar_cascade.py','follow_raar_circles.py','search_raar_circles.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    parts=[[] for _ in range(6)];begin=time.time()
    for offset in range(0,len(params),24):
        result=batch(params[offset:offset+24],100000,400000,20)
        np.savez_compressed(args.out/f'block-{offset:04}.npz',ids=np.arange(offset,min(offset+24,len(params))),params=params[offset:offset+24],rates=result[0],detlogs=result[1],bounds=result[2],ends=result[3],tails=result[4],counters=result[5])
        for a,r in zip(parts,result):a.append(r)
        print(json.dumps(dict(completed=min(offset+24,len(params)),total=len(params),seconds=time.time()-begin)),flush=True)
    rates,detlogs,bounds,ends,tails,counters=[np.concatenate(a) for a in parts];period=np.zeros(len(params),dtype=np.int64)
    for i,z in enumerate(tails):
        if counters[i,0]:continue
        for lag in range(1,513):
            if np.max(np.abs(z[lag:]-z[:-lag]))<1e-9:period[i]=lag;break
    flags=(counters[:,0]==0)&(rates[:,:,0].min(1)>.001)&(period==0)
    np.savez_compressed(args.out/'scan.npz',ids=np.arange(len(params)),params=params,rates=rates,detlogs=detlogs,bounds=bounds,ends=ends,tails=tails,counters=counters,period=period,flags=flags)
    report=dict(protocol,seconds=time.time()-begin,complete=int(np.sum(counters[:,0]==0)),positive=int(flags.sum()),period_counts={str(k):int(np.sum(period==k)) for k in np.unique(period)},max_trace_residual=float(np.max(np.abs(rates.sum(2)-detlogs))))
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

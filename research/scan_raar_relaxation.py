"""Declared relaxation slices through three selected RAAR circle regimes."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from follow_raar_circles import batch


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(8)
    with np.load(args.source) as f:source=f['params']
    params=[];groups=[]
    for ident in [1051,1483,820]:
        base=source[ident]
        for offset in np.linspace(-.02,.02,257):
            row=base.copy();row[0]+=offset;assert 0<row[0]<1
            params.append(row);groups.append(str(ident))
    params=np.array(params);groups=np.array(groups)
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),base_ids=[1051,1483,820],count=len(params),halfwidth=.02,points_per_slice=257,warm=100000,observed=400000,windows=20,initial='Fresh actual x=(.317,-.211,.419) for every setting',scope='Finite independent-start parameter slices, not a structural-stability proof or certified bifurcation diagram.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    for n in ['scan_raar_relaxation.py','follow_raar_circles.py','search_raar_circles.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    parts=[[] for _ in range(6)];begin=time.time()
    for offset in range(0,len(params),24):
        result=batch(params[offset:offset+24],100000,400000,20)
        np.savez_compressed(args.out/f'block-{offset:04}.npz',ids=np.arange(offset,min(offset+24,len(params))),params=params[offset:offset+24],groups=groups[offset:offset+24],rates=result[0],detlogs=result[1],bounds=result[2],ends=result[3],tails=result[4],counters=result[5])
        for a,r in zip(parts,result):a.append(r)
        print(json.dumps(dict(completed=min(offset+24,len(params)),total=len(params),seconds=time.time()-begin)),flush=True)
    rates,detlogs,bounds,ends,tails,counters=[np.concatenate(a) for a in parts];period=np.zeros(len(params),dtype=np.int64)
    for i,z in enumerate(tails):
        if counters[i,0]:continue
        for lag in range(1,513):
            if np.max(np.abs(z[lag:]-z[:-lag]))<1e-9:period[i]=lag;break
    flags=(counters[:,0]==0)&(rates[:,:,0].min(1)>.001)&(period==0)
    two=flags&(rates[:,:,1].min(1)>0)
    np.savez_compressed(args.out/'scan.npz',ids=np.arange(len(params)),params=params,groups=groups,rates=rates,detlogs=detlogs,bounds=bounds,ends=ends,tails=tails,counters=counters,period=period,flags=flags,two_positive=two)
    rows=[]
    for g in ['1051','1483','820']:
        m=groups==g;rows.append(dict(group=g,count=int(m.sum()),complete=int(np.sum(m&(counters[:,0]==0))),positive=int(np.sum(m&flags)),two_positive=int(np.sum(m&two)),periodic=int(np.sum(m&(period>0))),two_positive_beta_bounds=[float(params[m&two,0].min()),float(params[m&two,0].max())] if np.any(m&two) else None))
    report=dict(protocol,seconds=time.time()-begin,groups=rows,base_rows=[128,385,642],max_determinant_identity=float(np.max(np.abs(rates.sum(2)-detlogs))))
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

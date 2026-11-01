"""Continue every Yogi scout flag; preserve transient collapses and all tails."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from yogi_map import scout


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    d=np.load(args.source);a={k:d[k] for k in d.files}
    mask=(a['status']==0)&(a['windows'].min(1)>.005)&(np.ptp(a['tails'][:,:,0],axis=1)>1e-4)&(a['periods']==0)
    ids=np.flatnonzero(mask);p=a['parameters'][ids];starts=a['ends'][ids]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),
        warmup=0,observed=args.observed,window=50000,tail=4096,selection='Every original primary scout flag, without visual selection.',
        scope='Smooth-branch leading growth only. The map is discontinuous. No count of distinct chaotic attractors or priority assertion.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for n in ['follow_yogi.py','yogi_map.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    began=time.perf_counter();w,t,e,status,extrema,margins,changes=scout(p,starts,0,args.observed,50000,4096)
    periods=np.zeros(len(p),dtype=int)
    for i in range(len(p)):
        if status[i]!=0:continue
        tol=1e-9*(1+np.max(abs(t[i])))
        for lag in range(1,129):
            if np.max(abs(t[i,lag:]-t[i,:-lag]))<tol:periods[i]=lag;break
    np.savez_compressed(out/'histories.npz',ids=ids,parameters=p,starts=starts,windows=w,tails=t,ends=e,status=status,extrema=extrema,minimum_margins=margins,branch_changes=changes,periods=periods)
    keep=(status==0)&(w.min(1)>.005)&(np.ptp(t[:,:,0],axis=1)>1e-4)&(periods==0)
    records=[dict(id=int(ident),parameters=p[i].tolist(),minimum_window=float(w[i].min()),mean=float(w[i].mean()),
        last_window=float(w[i,-1]),tail_extent=np.ptp(t[i],axis=0).tolist(),period=int(periods[i]),retained=bool(keep[i])) for i,ident in enumerate(ids)]
    result=dict(declared=len(ids),completed=int(np.sum(status==0)),retained=int(keep.sum()),retained_ids=ids[keep].tolist(),
        near_zero_tails=int(np.sum(np.max(abs(t[:,:,:2]),axis=(1,2))<1e-100)),records=records,wall_seconds=time.perf_counter()-began)
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in result if k!='records'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--observed',type=int,default=1000000);main(p.parse_args())

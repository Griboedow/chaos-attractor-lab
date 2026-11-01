"""Continue every declared far endpoint and the weak full-volume check."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.spatial import cKDTree
from adamax_map import spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    source=np.load(args.source/'case-4099.npz');volume=np.load(args.source/'volume-checks.npz')
    far=np.flatnonzero(source['endpoint_distances']>.01);weak=np.flatnonzero(np.any(volume['windows'][:,:,0]<=.005,axis=1))
    starts=np.r_[source['snapshots'][-1,far],volume['ends'][weak]]
    parameters=np.r_[np.tile(source['parameters'],(len(far),1)),volume['parameters'][weak]]
    labels=[dict(kind='far_endpoint',index=int(i)) for i in far]+[dict(kind='weak_volume_check',index=int(i)) for i in weak]
    protocol=dict(source=str(args.source.resolve()),threshold=.01,far_count=len(far),weak_count=len(weak),warmup=50000,observed=500000,window=50000,tail=4096,
        selection='Every 2048-step endpoint farther than .01 normalized distance from the reference, plus every declared volume check with any leading window <=.005. No subsampling.',
        scope='Finite continuation of outliers. Regular closures and nonpositive windows remain negative outcomes.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['follow_adamax_outliers.py','adamax_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    w,tails,ends,resets,margins=spectra(parameters,starts,50000,500000,50000,4096,False)
    assert np.isfinite(w).all() and np.isfinite(tails).all()
    np.savez_compressed(out/'histories.npz',parameters=parameters,starts=starts,windows=w,tails=tails,ends=ends,reset_counts=resets,minimum_margins=margins)
    reference=np.load(args.reference,mmap_mode='r');tree=cKDTree(reference);scale=np.linalg.norm(np.ptp(reference,axis=0));rows=[]
    unique,inverse=np.unique(tails.reshape(-1,3),axis=0,return_inverse=True)
    distances=(tree.query(unique)[0]/scale)[inverse].reshape(tails.shape[:2])
    for i,label in enumerate(labels):
        tail=tails[i,-1024:];returns=[float(np.max(np.linalg.norm(tail[k:]-tail[:-k],axis=1))) for k in range(1,257)]
        period=next((k+1 for k,r in enumerate(returns) if r<1e-9*(1+np.linalg.norm(np.ptp(tail,axis=0)))),None)
        dist=distances[i]
        rows.append(dict(**label,spectrum_two=w[i].mean(0).tolist(),minimum_leading_window=float(w[i,:,0].min()),
            tail_period_up_to256=period,minimum_return_residual=min(returns),minimum_return_lag=int(np.argmin(returns)+1),
            support_distance_q50_q99_max=np.quantile(dist,[.5,.99,1]).tolist()))
    result=dict(runs=rows,completed=len(rows),declared=len(labels),positive_all_windows=sum(r['minimum_leading_window']>.005 for r in rows),
        closed_period_counts={str(p):sum(r['tail_period_up_to256']==p for r in rows) for p in sorted({r['tail_period_up_to256'] for r in rows if r['tail_period_up_to256'] is not None})},
        source_hashes={name:hashlib.sha256((args.source/name).read_bytes()).hexdigest() for name in ['case-4099.npz','volume-checks.npz']})
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print({k:result[k] for k in ['completed','positive_all_windows','closed_period_counts']},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

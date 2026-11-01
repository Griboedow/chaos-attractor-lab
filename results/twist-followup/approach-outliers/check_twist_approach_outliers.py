"""Continue every medusa ensemble endpoint farther than .02 from the reference."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.spatial import cKDTree
from twist_map import spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    d=np.load(args.source);ids=np.flatnonzero(d['nearest_reference_distance'][-1]>.02)
    p=np.tile(d['parameters'],(len(ids),1));starts=d['snapshots'][-1,ids]
    protocol=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        selection='Every candidate3254 endpoint whose Euclidean distance from archived reference tail is>.02 after4096iterates.',
        ids=ids.tolist(),count=len(ids),warmup=20000,observed=500000,window=50000,tail=20000,
        scope='Resolve unusual ensemble endpoints without relabeling their earlier4096step histories. No basin or asymptotic proof.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_twist_approach_outliers.py','twist_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    windows,tail,ends,means,radii=spectra(p,starts,20000,500000,50000,20000)
    path=out/'histories.npz';np.savez_compressed(path,ids=ids,parameters=p,starts=starts,windows=windows,tail=tail,ends=ends,means=means,radii=radii)
    ref=np.load(args.reference);i=list(ref['ids']).index(3254);tree=cKDTree(ref['tail'][i]);rows=[]
    for k,ident in enumerate(ids):
        distances=tree.query(tail[k],workers=4)[0]
        rows.append(dict(id=int(ident),spectrum=windows[k].mean(0).tolist(),minimum_leading_window=float(windows[k,:,0].min()),
            reference_distance_quantiles=np.quantile(distances,[.5,.95,.99,1]).tolist()))
    (out/'analysis.json').write_text(json.dumps(dict(count=len(ids),positive_all_windows=sum(r['minimum_leading_window']>0 for r in rows),
        maximum_trace_error=float(np.max(abs(windows.sum(2)-3*np.log(p[:,0,None])))),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),runs=rows),indent=2)+'\n')
    print('Completed',len(ids),'allpositive',sum(r['minimum_leading_window']>0 for r in rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

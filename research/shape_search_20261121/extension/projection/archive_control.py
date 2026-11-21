"""Apply matched consecutive-window descriptors to the frozen 108-regime archive."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from survey import BASE
from select_shapes import descriptors
ROOT=BASE.parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from data_cache import records,data_path


def windows(z):
    return np.stack([descriptors(z[n:n+1536]) for n in np.linspace(0,len(z)-1536,4,dtype=int)])


def main(args):
    prefix='results/attractor-novelty/reference/trajectories/'
    names=sorted({p.relative_to(ROOT).as_posix() for p in (ROOT/prefix).glob('*.npz')}|{p for p in records() if p.startswith(prefix)})
    assert len(names)==108
    with np.load(args.shortlist) as f:affine_spread=f['spread']
    with np.load(args.warped_shortlist) as f:warped_spread=f['spread']
    extracted=[];labels=[]
    for name in names:
        with np.load(data_path(name)) as f:z=f['points']
        x=windows(z)
        if np.isfinite(x).all():extracted.append(x);labels.append(Path(name).stem)
    ref=np.stack(extracted)
    out=[]
    for p in sorted(args.records.glob('case-*.json')):
        d=json.loads(p.read_text())
        if not d['retained']:continue
        with np.load(args.cache_root/d['cloud_relative']) as f:x=windows(f['points'])
        target=x.mean(0);m=ref.mean(1);spread=warped_spread if d['parameters'].get('warp',0) else affine_spread
        norms=lambda a:np.sqrt(np.mean((a/spread)**2,axis=-1))
        ds=norms(m-target);j=int(np.argmin(ds));own=float(norms(x-target).max())
        archive_spread=float(norms(ref[j]-m[j]).max())
        out.append(dict(case=d['case'],nearest=labels[j],distance=float(ds[j]),self_max_distance=own,nearest_reference_max_distance=archive_spread))
    args.out.write_text(json.dumps(dict(archive_total=len(names),finite=len(ref),method='Four matched contiguous 1536-state windows and the same feature scaling as the selection. PCA descriptors are numerical morphology indicators.',rows=out),indent=2)+'\n',encoding='utf8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--records',type=Path,required=True);p.add_argument('--cache-root',type=Path,required=True);p.add_argument('--shortlist',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--warped-shortlist',type=Path,required=True)
    main(p.parse_args())

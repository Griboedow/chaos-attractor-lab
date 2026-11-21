"""Compare finite persistence diagrams with the 108 frozen catalog regimes."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
import gudhi
from topology import analyze,save

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools'))
from data_cache import data_path,records


def distance(a,b):
    value=lower=upper=0.
    for k in ('0','1','2'):
        x=np.array(a[k]['bars'],float).reshape(-1,2)
        y=np.array(b[k]['bars'],float).reshape(-1,2)
        d=float(gudhi.bottleneck_distance(x,y,e=1e-9))
        ex=.5*(x[-1,1]-x[-1,0]) if len(x)==6 else 0.
        ey=.5*(y[-1,1]-y[-1,0]) if len(y)==6 else 0.
        value=max(value,d)
        lower=max(lower,d-ex-ey-1e-9)
        upper=max(upper,d+ex+ey+1e-9)
    return np.array([value,lower,upper])


def main(args):
    args.cache.mkdir(parents=True,exist_ok=True)
    base='results/attractor-novelty/reference/trajectories/'
    files={p.relative_to(ROOT).as_posix() for p in (ROOT/base).glob('*.npz')}
    files.update(k for k in records() if k.startswith(base))
    assert len(files)==108
    refs=[];omitted=[]
    for j,name in enumerate(sorted(files)):
        key=Path(name).stem
        target=args.cache/f'{key}.json'
        if not target.exists():
            source=data_path(name)
            with np.load(source) as f:z=f['points'].astype(float)
            scale=float(np.sqrt(np.sum((z-z.mean(0))**2,axis=1).mean()))
            if scale<1e-10:
                save(target,dict(case=key,skipped='RMS radius below 1e-10',rms_radius=scale))
            else:
                row=analyze(z,key,9134+j,[1536,3072])
                row['sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
                save(target,row)
        row=json.loads(target.read_text())
        if 'skipped' in row: omitted.append(row)
        else: refs.append(row)
    results=[]
    for path in sorted(args.candidates.glob('persistence-*.json')):
        row=json.loads(path.read_text())
        sample=row['replicates']
        matches=[]
        for ref in refs:
            scores=[distance(a['groups'],b['groups']) for a,b in zip(sample,ref['replicates'])]
            value,lo,hi=np.median(scores,axis=0)
            matches.append(dict(reference=ref['case'],truncated_median=float(value),full_median_lower=float(lo),full_median_upper=float(hi)))
        matches.sort(key=lambda x:x['truncated_median'])
        # Matched sample sizes in opposite temporal halves.
        self_scores=[distance(sample[j]['groups'],sample[j+2]['groups']).tolist() for j in (0,1)]
        results.append(dict(case=row['case'],matches=matches[:4],self_by_size=self_scores,smallest_full_median_lower_across_refs=min(m['full_median_lower'] for m in matches)))
    study=dict(reference_count=len(refs),skipped=omitted,gudhi=gudhi.__version__,method='Within each of four matched finite samples, compute the maximum bottleneck distance in H0, H1 and H2 with the six longest finite radius intervals per degree. Omitted bars can be paired to their diagonals within half the sixth lifetime. Triangle bounds propagate this truncation to the full finite diagrams. Rank by median over four paired samples. These are point-cloud offsets in RMS units.',self_columns=['truncated','full_lower','full_upper'],cases=results)
    save(args.out,study)
    print(len(refs),'reference clouds;',len(results),'candidate clouds',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--candidates',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())

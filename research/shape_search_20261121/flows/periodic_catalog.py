"""Index distinct accepted stable-cycle cases and check for repeated traversals."""
import argparse
import csv
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np


def shifts(points):
    n=len(points)
    scale=np.linalg.norm(np.std(points,axis=0))
    phases=(np.linspace(0,n,19,endpoint=False)+n*.013)%n
    def at(x):
        k=np.floor(x).astype(int)
        f=(x-k)[:,None]
        return points[k%n]*(1-f)+points[(k+1)%n]*f
    start=at(phases)
    # Any traversal repeated k times would agree after exactly one kth of its period at every phase.
    residual=[np.max(np.linalg.norm(at(phases+n/k)-start,axis=1))/scale for k in range(2,257)]
    return float(min(residual)),int(np.argmin(residual))+2


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    manifests=json.loads((a.out/'manifest.json').read_text())
    selected={r['case'] for r in manifests['stable_cycles']}
    files={}
    # Full screening takes precedence over the earlier overlapping pilot.
    for folder in sorted(a.cache.glob('knot-validated-*-*'),key=lambda p:('pre' not in p.name,p.name)):
        for file in sorted(folder.glob('case-*.json')):
            files[int(file.stem.split('-')[1])]=file
    result=[];records=[];signatures={}
    for case,file in sorted(files.items()):
        info=json.loads(file.read_text());cloud=file.with_suffix('.npz')
        with np.load(cloud) as data:
            pts=data['points'];floor,k=shifts(pts)
        if floor<.0001:raise ValueError(f'Possible multiple traversal of case {case} at T/{k}')
        coeff=info['alexander_coefficients'];key=tuple(coeff)
        signatures.setdefault(key,[]).append(case)
        record=dict(case=case,family=info['family'],selected=case in selected,period=info['period'],
                    dominant_transverse_modulus=info['max_transverse_multiplier_abs'],
                    closure_rms=info['independent_closure_rms_units'],determinant=info['determinant'],
                    alexander_coefficients=json.dumps(coeff,separators=(',',':')),
                    min_crossings_seen=min(v['crossings'] for v in info['views']),
                    returning_1pct=sum(t['p95_distance']<.02 for t in info['recovery'][:3]),
                    returning_4pct=sum(t['p95_distance']<.02 for t in info['recovery'][3:]),
                    minimum_phase_shift_rms=floor,worst_divisor=k,equation_seed=info['equation_seed'],
                    refined_start=json.dumps(info['refined_start'],separators=(',',':')),points=len(pts),
                    sha256=hashlib.sha256(cloud.read_bytes()).hexdigest(),
                    cloud_relative=cloud.relative_to(a.cache.parent).as_posix())
        result.append(record)
        records.append(dict(case=case,family=info['family'],collection='stable_periodic_orbit',
            label=f"Periodic knot, determinant {info['determinant']}",path=str(cloud.resolve()),
            cloud_relative=record['cloud_relative'],sha256=record['sha256'],point_key='points',
            shape=[len(pts),3],point_dtype='float64',equation_seed=info['equation_seed'],
            parameters=info['parameters'],initial=info['refined_start'],period=info['period'],
            sample_time_step=info['period']/len(pts),alexander_coefficients=coeff,
            dominant_transverse_modulus=info['max_transverse_multiplier_abs']))
    with (a.out/'periodic-catalog.csv').open('w',newline='',encoding='utf-8') as handle:
        w=csv.DictWriter(handle,fieldnames=result[0]);w.writeheader();w.writerows(result)
    (a.out/'periodic-clouds.jsonl').write_text(''.join(json.dumps(row,separators=(',',':'))+'\n' for row in records))
    screening={}
    for folder in sorted(a.cache.glob('knots-*')):
        with (folder/'screen.csv').open() as f:
            counter=Counter(r['status'] for r in csv.DictReader(f))
        screening[folder.name]=dict(count=sum(counter.values()),outcomes=counter,
                                   run=json.loads((folder/'run.json').read_text()))
    attempts={}
    for folder in sorted(a.cache.glob('knot-validated-*-*')):
        path=folder/'attempts.json'
        if path.exists():attempts[folder.name]=dict(Counter(json.loads(path.read_text()).values()))
    summary=dict(screening=screening,validation_attempts=attempts,distinct_accepted_parameter_cases=len(result),
                 families=Counter(r['family'] for r in result),different_sampled_Alexander_polynomials=len(signatures),
                 primitive_control=dict(phases=19,divisors=[2,256],minimum_of_all_phase_shift_rms=min(r['minimum_phase_shift_rms'] for r in result)),
                 signatures=[dict(coefficients=list(key),cases=val) for key,val in sorted(signatures.items(),key=lambda group:(len(group[0]),group[0]))])
    (a.out/'periodic-accounting.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('cases',len(result),'polynomials',len(signatures),'families',dict(summary['families']))


if __name__=='__main__':
    main()

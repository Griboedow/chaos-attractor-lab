"""Package all completed doubled orbits at the exact final chaotic coefficients."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from cycle_scan import resample
from exact_knot import prepare
from polygon_isotopy import calculate
from symmetry_links import link_integer


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--packet',type=Path,required=True)
    ap.add_argument('--cache',type=Path,required=True)
    a=ap.parse_args();r=a.packet
    m=json.loads((r/'manifest.json').read_text())
    chaos=next(row for row in m['cases'] if row['case']==180844)
    main=next(v for v in m['knotted_periodic_skeletons'] if v['case']==180844)
    p=np.array(chaos['parameters'])
    reports=[];stored=[];low={1:a.cache/main['cloud_relative']}
    for file in sorted(a.source.glob('order-*.json'),key=lambda f:int(f.stem[6:])):
        d=json.loads(file.read_text());reports.append(d)
        if not d['completed']:continue
        path=a.cache/d['cloud_relative'];n=d['order']
        child=np.load(path)['points'];inf=d['final']
        record=dict(case=180844,period_order=n,family='ring',label=f'Unstable periodic orbit of order {n}',
                    parameters=p.tolist(),initial=d['refined_start'],period=d['period'],
                    dominant_multiplier=inf['dominant'],neutral_multiplier=inf['neutral'],
                    closure_native=inf['closure'],orbit_to_cloud_p95_rms=d['cloud_distance_p95'],
                    cloud_to_orbit_p95_rms=d['reverse_distance_p95'],shape=list(child.shape),
                    point_key='points',point_dtype='float64',sample_time_step=d['period']/len(child),
                    path=str(path.resolve()),cloud_relative=d['cloud_relative'],sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        if n==2:
            expected=json.loads((r/'doubling-cable.json').read_text())['child_alexander']
            witness=calculate(path,expected,True);witness['source_relative']=d['cloud_relative']
            (r/'nested-dense-witness.json').write_text(json.dumps(witness,separators=(',',':'))+'\n')
            record['dense_polygon_witness']='nested-dense-witness.json'
            record['alexander_coefficients']=expected
        if n//2 in low:
            parent=np.load(low[n//2])['points'];scale=np.linalg.norm(np.std(parent,axis=0));links=[]
            for size in [2560,5120,16384,len(child)]:
                first=resample(parent,size) if size!=len(child) else parent
                second=resample(child,size) if size!=len(child) else child
                pair=[np.rint(v*1e9/scale).astype(np.int64) for v in (first,second)]
                checks=[]
                for seed in (17,35,81):
                    _,mat,_,_=prepare(first[::32],seed);value,cross=link_integer(pair[0],pair[1],mat)
                    checks.append([value,cross])
                links.append(dict(parent_vertices=len(first),child_vertices=len(second),values_and_crossings=checks))
            record['adjacent_core_pair']=dict(parent_order=n//2,camera_seeds=[17,35,81],link_checks=links)
        low[n]=path;stored.append(record)
    (r/'nested-skeletons.json').write_text(json.dumps(dict(case=180844,cores=stored),indent=2)+'\n')
    (r/'nested-continuation.json.gz').write_bytes(gzip.compress(json.dumps(reports,separators=(',',':')).encode(),mtime=0))
    print([(v['period_order'],v['dominant_multiplier'],v.get('adjacent_core_pair',{}).get('link_checks',[{}])[-1]) for v in stored])


if __name__=='__main__':
    main()

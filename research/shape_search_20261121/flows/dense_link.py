"""Compare exact linking under mesh refinement for closely spaced period doublings."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from cycle_scan import resample
from exact_knot import prepare
from symmetry_links import link_integer,check_controls


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pairs',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True);check_controls()
    rows=json.loads((a.pairs/'report.json').read_text());log=[]
    for record in rows:
        k=record['child_period_order']
        both=np.load(a.pairs/f'pair-{k}.npz')['points'];core,child=np.array_split(both,2)
        scale=np.linalg.norm(np.std(core,axis=0))
        info=dict(order=k,scale=scale,meshes=[])
        meshes=[2560,5120,10240,20480,32768]
        if len(core)>32768:meshes.extend([49152,len(core)])
        for n in meshes:
            for spacing in (['arclength','time'] if n==len(core) else ['arclength']):
                closed=[(resample(v,n) if spacing=='arclength' else v) for v in (core,child)]
                polys=[np.rint(v/scale*1e9).astype(np.int64) for v in closed]
                views=[]
                for seed in (17,35,81):
                    _,m,_,_=prepare(core[::64],seed)
                    num,cross=link_integer(polys[0],polys[1],m)
                    views.append(dict(camera=m.tolist(),linking=num,crossings=cross))
                assert len({view['linking'] for view in views})==1
                info['meshes'].append(dict(vertices=n,spacing=spacing,views=views))
                print(k,n,spacing,views[0]['linking'],flush=True)
            (a.out/f'order-{k}.json').write_text(json.dumps(info,indent=2)+'\n')
        log.append(info)
        (a.out/'report.json').write_text(json.dumps(log,indent=2)+'\n')


if __name__=='__main__':
    main()

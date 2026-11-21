"""Test which symmetry-related attracting cycle captures each perturbed start."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.spatial import cKDTree
from field import integrate,velocity_tangent,ring_tangent


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--links',type=Path,required=True)
    ap.add_argument('--source',type=Path,nargs='+',required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    for cert in args.links.glob('case-*.json'):
        w=json.loads(cert.read_text());case=w['case']
        record=next(json.loads((d/f'case-{case}.json').read_text()) for d in args.source if (d/f'case-{case}.json').exists())
        p=np.asarray(record['parameters']);period=record['period'];rms=w['common_rms']
        rhs=velocity_tangent if record['family']=='cartesian' else ring_tangent
        with np.load(cert.with_suffix('.npz')) as f:
            all_points=f['points'];offsets=f['offsets']
        clouds=[all_points[a:b] for a,b in zip(offsets,offsets[1:])]
        trees=[cKDTree(v) for v in clouds]
        starts=[]
        rng=np.random.default_rng(493189+case)
        npoints=len(clouds[0])
        for i in range(6):
            target=clouds[0][i*npoints//6]
            noise=rng.normal(size=3);radius=.01 if i<3 else .04
            starts.append((0,target+noise/np.linalg.norm(noise)*radius*rms,radius))
        rng=np.random.default_rng(28811+case)
        for j,v in enumerate(clouds):
            for i in range(3):
                target=v[(i*8719+921*j)%npoints]
                noise=rng.normal(size=3);radius=.005
                starts.append((j,target+noise/np.linalg.norm(noise)*radius*rms,radius))
        rows=[]
        for source,start,radius in starts:
            h=min(.0125,period/800)
            _,tail=integrate(p,start,h,round(24*period/h),round(3*period/h),1,max(1,round(period/250/h)),rhs)
            distances=np.array([np.quantile(tree.query(tail,workers=1)[0],.95)/rms for tree in trees])
            winner=int(np.argmin(distances))
            rows.append(dict(source=source,initial_xyz=start.tolist(),radius_rms=radius,distances=distances.tolist(),
                             nearest_component=winner,captured=bool(distances[winner]<.02)))
        equiv=0.
        for mat in map(np.asarray,w['transforms']):
            for v in clouds[0][::1321]:
                z=np.r_[v,0.,0.,0.]
                transformed=np.r_[v@mat,0.,0.,0.]
                equiv=max(equiv,float(np.max(np.abs(rhs(transformed,p)[:3]-rhs(z,p)[:3]@mat))))
        out=dict(case=case,component_count=len(clouds),trials=rows,maximum_equivariance_discrepancy=equiv,
                 independent_fixed_step='RK4, 24 settling periods and 3 recorded periods')
        (args.out/cert.name).write_text(json.dumps(out,indent=2)+'\n')
        print(case,[r['nearest_component'] if r['captured'] else None for r in rows],flush=True)


if __name__=='__main__':
    main()

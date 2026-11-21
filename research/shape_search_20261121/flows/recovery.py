"""Displace states along a late trajectory and test finite-time recovery."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.spatial import cKDTree
from field import integrate,velocity_tangent,ring_tangent,rotor_tangent,plane_tangent,crossed_tangent


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True,type=Path)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--ids',type=int,nargs='+',required=True)
    ap.add_argument('--family',choices=('cartesian','ring','rotor','planes','crossed'),required=True)
    ap.add_argument('--radii',type=float,nargs='+',default=[.01,.05])
    ap.add_argument('--starts',type=int,default=8)
    ap.add_argument('--dt',type=float,default=.015)
    ap.add_argument('--burn',type=float,default=200.)
    ap.add_argument('--window',type=float,default=250.)
    ap.add_argument('--blocks',type=int,default=6)
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    (args.out/'run.json').write_text(json.dumps({**vars(args),'out':str(args.out),'source':str(args.source)},indent=2)+'\n')
    rhs=dict(cartesian=velocity_tangent,ring=ring_tangent,rotor=rotor_tangent,planes=plane_tangent,crossed=crossed_tangent)[args.family]
    for case in args.ids:
        with np.load(args.source/f'case-{case}.npz') as f:
            z=f['points'];p=f['parameters']
        scale=np.sqrt(np.mean(np.sum((z-z.mean(0))**2,axis=1)))
        tree=cKDTree(z)
        rng=np.random.default_rng(np.random.SeedSequence([20261121,case,425]))
        records=[]
        for radius in args.radii:
            group=[]
            for i in range(args.starts):
                loc=int((i+rng.uniform(.15,.85))*len(z)/args.starts)
                vector=rng.normal(size=3)
                start=z[loc]+vector/np.linalg.norm(vector)*radius*scale
                rates,tail=integrate(p,start,args.dt,round(args.burn/args.dt),round(args.window/args.dt),
                                     args.blocks,max(1,round(.5/args.dt)),rhs)
                distance=tree.query(tail[::2],workers=1)[0]/scale
                q=float(np.quantile(distance,.95))
                group.append(dict(index=loc,initial_xyz=start.tolist(),rates=rates.tolist(),
                                  distance_p95_rms_units=q,pass_growth_and_distance=bool(min(rates)>.008 and q<.06)))
            records.append(dict(radius_rms=float(radius),starts=len(group),passed=sum(v['pass_growth_and_distance'] for v in group),
                                positive_windows=sum(min(v['rates'])>.008 for v in group),trials=group))
        (args.out/f'case-{case}.json').write_text(json.dumps(dict(case=case,family=args.family,groups=records),indent=2)+'\n')
        print(case,[v['passed'] for v in records],flush=True)


if __name__=='__main__':
    main()

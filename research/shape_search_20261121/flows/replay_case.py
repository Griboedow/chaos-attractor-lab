"""Regenerate one final native XYZ cloud solely from the portable case record."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.integrate import solve_ivp
from adaptive import full_rhs
from field import integrate,velocity_tangent,ring_tangent,rotor_tangent,plane_tangent,crossed_tangent
from field import make_case,make_ring,make_rotor,make_planes,make_crossed


def orbit(start,period,p,mode,n,tolerance):
    sol=solve_ivp(full_rhs,(0.,period),np.r_[start,np.eye(3).ravel(),0.],args=(p,mode),
                  method='DOP853',rtol=tolerance,atol=tolerance/100,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return sol.sol(np.linspace(0,period,n,endpoint=False))[:3].T


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',type=Path,default=Path('results/shape_search_20261121/flows'))
    ap.add_argument('--case',type=int,required=True)
    ap.add_argument('--out',type=Path,required=True)
    group=ap.add_mutually_exclusive_group()
    group.add_argument('--linked',action='store_true',help='Generate all symmetry-related members with component labels')
    group.add_argument('--skeleton',action='store_true',help='Generate the unstable knotted orbit of a chaotic band')
    group.add_argument('--doubling',type=int,choices=(2,4,8,16),help='Generate one of the ring period-doubling pairs')
    group.add_argument('--skeleton-order',type=int,help='Generate a higher-period unstable core at the final chaotic coefficients')
    args=ap.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.parent.mkdir(exist_ok=True,parents=True)
    manifest=json.loads((args.packet/'manifest.json').read_text())
    stable={row['case']:row for row in manifest.get('stable_cycles',[])}
    if not args.skeleton and args.case not in stable and (args.packet/'periodic-catalog.csv').exists():
        with (args.packet/'periodic-catalog.csv').open() as stream:
            archive=next((r for r in csv.DictReader(stream) if int(r['case'])==args.case),None)
        if archive is not None:
            generator=dict(cartesian=make_case,ring=make_ring,rotor=make_rotor,planes=make_planes,crossed=make_crossed)[archive['family']]
            p,_=generator(int(archive['equation_seed']),args.case)
            stable[args.case]=dict(family=archive['family'],parameters=p.tolist(),
                initial=json.loads(archive['refined_start']),period=float(archive['period']),
                shape=[int(archive['points']),3],sha256=archive['sha256'])
    if args.linked and not any(r['case']==args.case for r in manifest['linked_cycles']):
        ap.error('No linked collection is registered for this case')
    if args.skeleton_order:
        rows=json.loads((args.packet/'nested-skeletons.json').read_text())['cores']
        row=next((x for x in rows if x['case']==args.case and x['period_order']==args.skeleton_order),None)
        if row is None:ap.error('No matching continued orbit at the final coefficients')
        if row.get('method')=='C7 relative shooting':
            from relative_cores import sector_cloud
            points,_=sector_cloud(np.array(row['initial']),row['sector_period'],np.array(row['parameters']),
                                  row['relative_rotation_index'],row['shape'][0])
        else:points=orbit(row['initial'],row['period'],np.array(row['parameters']),1,row['shape'][0],3e-12)
        np.savez_compressed(args.out,points=points)
    elif args.doubling:
        if args.case!=180844:ap.error('Doubled-orbit pairs belong to the Case 180844 parameter path')
        row=next(r for r in json.loads((args.packet/'doubling-pairs.json').read_text()) if r['child_period_order']==args.doubling)
        p=np.array(row['parameters']);n=row['shape'][0]//2
        clouds=[orbit(np.array(row[k+'_start']),row[k]['period'],p,1,n,3e-12) for k in ('parent','child')]
        points=np.concatenate(clouds)
        np.savez_compressed(args.out,points=points,component=np.repeat([0,1],n),offsets=[0,n,2*n])
    elif args.case in stable or args.skeleton:
        if args.skeleton:
            row=next(r for r in manifest['knotted_periodic_skeletons'] if r['case']==args.case)
            parent=next(r for r in manifest['cases'] if r['case']==args.case)
            p=np.array(parent['parameters']).ravel();start=np.array(row['refined_start'])
        else:
            row=stable[args.case];p=np.array(row['parameters']).ravel();start=np.array(row['initial'])
        t=row['period'];mode=['cartesian','ring','rotor','planes','crossed'].index(row['family'])
        tolerance=3e-12 if args.skeleton else 2e-12
        points=orbit(start,t,p,mode,row['shape'][0],tolerance)
        if args.linked:
            row=next(r for r in manifest['linked_cycles'] if r['case']==args.case)
            clouds=[points@np.array(m) for m in row['row_vector_transforms']]
            offsets=np.arange(len(clouds)+1)*len(points)
            component=np.repeat(np.arange(len(clouds)),len(points))
            points=np.concatenate(clouds)
            np.savez_compressed(args.out,points=points,component=component,offsets=offsets)
        else:
            np.savez_compressed(args.out,points=points)
    else:
        row=next(r for r in manifest['cases'] if r['case']==args.case)
        spec=json.loads((args.packet/'run-settings.json').read_text())[row['dense_run']]
        p=np.array(row['parameters']);start=np.array(row['initial'])
        rhs={'cartesian':velocity_tangent,'ring':ring_tangent,'rotor':rotor_tangent,
             'planes':plane_tangent,'crossed':crossed_tangent}[row['family']]
        step=spec['dt']
        rates,points=integrate(p,start,step,round(spec['burn']/step),round(spec['window']/step),
                               spec['blocks'],spec['stride'],rhs)
        np.savez_compressed(args.out,points=points,parameters=p,initial=start,rates=rates,
                            step_size=step,sample_stride=spec['stride'])
    print(json.dumps(dict(case=args.case,shape=points.shape,sha256=hashlib.sha256(args.out.read_bytes()).hexdigest(),
                     original_sha256=row['sha256'])))


if __name__=='__main__':
    main()

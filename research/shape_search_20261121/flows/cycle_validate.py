"""Check shooting closure, attracting Floquet multipliers and polygon knot invariants."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
from adaptive import full_rhs
from cycle_scan import resample
from field import integrate,make_case,make_ring,make_rotor,make_planes,make_crossed
from field import velocity_tangent,ring_tangent,rotor_tangent,plane_tangent,crossed_tangent
from knot_diagram import alexander,positive_controls
from orbits import refine


def validate(hit,mode,gen,rhs,seed):
    p,generated=gen(seed,hit['case'])
    out=refine(hit['axis'],np.asarray(hit['start']),hit['period'],hit['rms'],p.ravel(),mode)
    if out is None:
        return None,'Newton iteration'
    start,period=out
    y=np.r_[start,np.eye(3).ravel(),0.]
    sol=solve_ivp(full_rhs,(0.,period),y,args=(p.ravel(),mode),method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
    if not sol.success:
        return None,'adaptive solver'
    end=sol.y[:,-1]
    closure=np.linalg.norm(end[:3]-start)/hit['rms']
    if closure>5e-7:
        return None,'independent closure'
    if min(np.linalg.norm(sol.sol(period/k)[:3]-start)/hit['rms'] for k in range(2,17))<5e-5:
        return None,'multiple traverse'
    multipliers=np.linalg.eigvals(end[3:12].reshape(3,3))
    k=np.argmin(abs(multipliers-1))
    rest=np.delete(multipliers,k)
    if abs(multipliers[k]-1)>1e-4 or max(abs(rest))>=.995:
        return None,'weak or unstable transverse spectrum'
    points=sol.sol(np.linspace(0,period,32768,endpoint=False))[:3].T
    views=[]
    for count in (320,640,1280):
        polygon=resample(points,count)
        group=[]
        for camera in (17,35,81,192,261,713,967,1217,2114,3117,5032,7101):
            diagram=alexander(polygon,camera)
            if diagram is not None:
                group.append((diagram['crossings'],camera,diagram))
        group.sort(key=lambda a:a[:2])
        if len(group)<3 or any(x[2]['determinant']!=hit['determinant'] for x in group):
            return None,'unstable diagram determinant'
        for _,camera,diagram in group[:2]:
            if diagram['crossings']>85:
                return None,'projection too complicated'
            polynomial=alexander(polygon,camera,True)
            views.append(dict(vertices=count,camera_seed=camera,**polynomial))
    coefficients=views[0]['coefficients']
    if any(v['coefficients']!=coefficients for v in views):
        return None,'unstable Alexander polynomial'
    tree=cKDTree(points)
    scale=np.sqrt(np.mean(np.sum((points-points.mean(0))**2,axis=1)))
    rng=np.random.default_rng(493189+hit['case'])
    trials=[]
    for i in range(6):
        idx=(i*len(points))//6
        noise=rng.normal(size=3)
        target_radius=(.01 if i<3 else .04)
        perturbed=points[idx]+noise/np.linalg.norm(noise)*target_radius*scale
        h=min(.0125,period/800)
        _,tail=integrate(p,perturbed,h,round(24*period/h),round(3*period/h),1,max(1,round(period/250/h)),rhs)
        distance=tree.query(tail,workers=1)[0]/scale
        trials.append(dict(relative_radius=target_radius,p95_distance=float(np.quantile(distance,.95))))
    record=dict(case=hit['case'],family=hit['family'],equation_seed=seed,parameters=p.tolist(),
                original_start=generated.tolist(),refined_start=start.tolist(),period=period,
                independent_closure_rms_units=float(closure),neutral_multiplier=float(multipliers[k].real),
                transverse_multipliers=[[float(x.real),float(x.imag)] for x in rest],
                max_transverse_multiplier_abs=float(max(abs(rest))),determinant=hit['determinant'],
                alexander_coefficients=coefficients,views=views,recovery=trials,
                winding_about_z=float((np.unwrap(np.arctan2(points[:,1],points[:,0]))[-1]-np.arctan2(points[0,1],points[0,0]))/2/np.pi))
    return (record,points),'accepted'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path)
    ap.add_argument('--hits',type=Path)
    ap.add_argument('--ids',type=int,nargs='+')
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--limit',type=int,default=80)
    ap.add_argument('--per-det',type=int,default=8)
    ap.add_argument('--resume',action='store_true',help='Skip previously recorded attempts in the output directory')
    ap.add_argument('--seed',type=int,default=2026112123)
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=args.resume);positive_controls()
    hit_path=args.hits or (args.source/'hits.jsonl')
    hits=[json.loads(s) for s in hit_path.read_text().splitlines() if s.strip()]
    if args.ids:
        hits=[h for h in hits if h['case'] in args.ids]
    # Prefer simpler witness projections and sample each determinant before its runners-up.
    by_det={}
    for hit in hits:
        by_det.setdefault(hit['determinant'],[]).append(hit)
    items=[]
    for rows in by_det.values():
        rows.sort(key=lambda r:(min(r['two_crossing_counts']),r['period'],r['case']))
    for i in range(args.per_det):
        for det,rows in sorted(by_det.items()):
            if i<len(rows):
                items.append(rows[i])
    items=items[:args.limit]
    report_file=args.out/'attempts.json'
    report=json.loads(report_file.read_text()) if args.resume and report_file.exists() else {}
    for hit in items:
        if str(hit['case']) in report:
            continue
        mode=('cartesian','ring','rotor','planes','crossed').index(hit['family'])
        gen=(make_case,make_ring,make_rotor,make_planes,make_crossed)[mode]
        rhs=(velocity_tangent,ring_tangent,rotor_tangent,plane_tangent,crossed_tangent)[mode]
        result,status=validate(hit,mode,gen,rhs,args.seed)
        report[str(hit['case'])]=status
        if result:
            record,points=result
            np.savez_compressed(args.out/f"case-{hit['case']}.npz",points=points)
            (args.out/f"case-{hit['case']}.json").write_text(json.dumps(record,indent=2)+'\n')
            print(hit['case'],record['alexander_coefficients'],flush=True)
        report_file.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()

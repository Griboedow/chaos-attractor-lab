"""Follow a stable knotted cycle through parameter changes and test its final orbit."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
from adaptive import full_rhs,polar_ring
from cycle_scan import resample
from exact_knot import exact_diagram,prepare
from orbits import refine


def dynamics(start, period, p, dense=False):
    sol=solve_ivp(full_rhs,(0.,period),np.r_[start,np.eye(3).ravel(),0.],args=(p,1),method='DOP853',
                  rtol=3e-12 if dense else 2e-10,atol=3e-14 if dense else 2e-12,dense_output=dense)
    if not sol.success:return None
    final=sol.y[:,-1]
    eig=np.linalg.eigvals(final[3:12].reshape(3,3));i=np.argmin(abs(eig-1))
    nontrivial=np.delete(eig,i)
    dominant=nontrivial[np.argmax(abs(nontrivial))]
    row=dict(period=float(period),neutral=[float(eig[i].real),float(eig[i].imag)],
             dominant=[float(dominant.real),float(dominant.imag)],
             modulus=float(abs(dominant)),closure=float(np.linalg.norm(final[:3]-start)))
    return row,sol


def continuation(original, p_final, start, period, scale):
    current=start.copy();t=period;position=0.;step=1/16
    initial=dynamics(current,t,original)
    rows=[dict(fraction=0.,**initial[0])]
    while position<1-1e-12:
        target=min(1.,position+step)
        p=(1-target)*original+target*p_final
        axis=int(np.argmax(abs(polar_ring(current,p))))
        got=refine(axis,current,t,scale,p,1)
        if got is None:
            step/=2
            if step<1/512:return rows,None
            continue
        candidate,new_t=got
        valid=dynamics(candidate,new_t,p)
        if valid is None or abs(valid[0]['neutral'][0]-1)>.01 or valid[0]['closure']/scale>1e-5:
            step/=2
            if step<1/512:return rows,None
            continue
        current=candidate;t=new_t
        rows.append(dict(fraction=target,**valid[0]))
        position=target
        step=min(1/12,1.2*step)
    return rows,(current,t)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--leads',type=Path,required=True)
    ap.add_argument('--bases',type=Path,required=True)
    ap.add_argument('--dense',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--ids',type=int,nargs='+',required=True)
    a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    leads={row['case']:row for row in json.loads(a.leads.read_text())}
    for case in a.ids:
        meta=leads[case];base=meta['base']
        record=json.loads((a.bases/f'case-{base}.json').read_text())
        p0=np.array(record['parameters']);start=np.array(record['refined_start']);t=record['period']
        with np.load(a.dense/f'case-{case}.npz') as f:
            p=f['parameters'];chaos=f['points']
        scale=np.linalg.norm(np.std(chaos,axis=0))
        rows,got=continuation(p0,p,start,t,scale)
        report=dict(case=case,base=base,branch=rows,completed=got is not None)
        if got is not None:
            new_start,new_t=got
            summary,sol=dynamics(new_start,new_t,p,True)
            orbit=sol.sol(np.linspace(0,new_t,32768,endpoint=False))[:3].T
            inner=cKDTree(chaos).query(orbit[::4])[0]/scale
            outer=cKDTree(orbit).query(chaos[::16])[0]/scale
            polygon=resample(orbit,640);views=[]
            for camera in [v['camera_seed'] for v in record['views'][:2]]:
                vertices,matrix,shift,units=prepare(polygon,camera)
                exact=exact_diagram(vertices,matrix);exact.pop('minimum_depth',None)
                views.append(dict(matrix=matrix.tolist(),**exact))
            report.update(refined_start=new_start.tolist(),final=summary,final_exact_views=views,
                          same_alexander_polynomial=all(v['coefficients']==record['alexander_coefficients'] for v in views),
                          orbit_to_cloud_distance_p95=float(np.quantile(inner,.95)),
                          cloud_to_orbit_distance_p95=float(np.quantile(outer,.95)),
                          projection_witness=dict(integer_vertices=vertices.tolist(),center=shift.tolist(),rms=units))
            np.savez_compressed(a.out/f'case-{case}.npz',points=orbit)
        (a.out/f'case-{case}.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
        print(case,report['completed'],report.get('final',{}).get('dominant'),flush=True)


if __name__=='__main__':
    main()

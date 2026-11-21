"""Replay the first observed period doublings on the fixed seven-crossing case path."""
import argparse
import hashlib
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.spatial import cKDTree
from band_continuation import dynamics
from cycle_scan import resample
from exact_knot import prepare
from orbits import refine
from symmetry_links import link_integer


def track(start,period,f0,f1,parameters,scale):
    times=np.linspace(f0,f1,max(2,int(np.ceil(abs(f1-f0)/.001))+1))[1:]
    for f in times:
        p=parameters(f)
        # Trust the previous seed's phase plane. Vary the free coordinates and time only.
        from adaptive import polar_ring
        axis=int(np.argmax(abs(polar_ring(start,p))))
        got=refine(axis,start,period,scale,p,1)
        if got is None:return None
        start,period=got
    return start,period


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',type=Path,default=Path('results/shape_search_20261121/flows'))
    ap.add_argument('--primary',type=Path,required=True)
    ap.add_argument('--long-return',type=Path)
    ap.add_argument('--orders',nargs='+',type=int,default=[2,4,8])
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    m=json.loads((a.packet/'manifest.json').read_text())
    all_rows=json.loads(a.primary.read_text())['trials'];rows={round(r['fraction'],4):r for r in all_rows}
    p0=np.array(next(r for r in m['stable_cycles'] if r['case']==28756)['parameters'])
    p1=np.array(next(r for r in m['cases'] if r['case']==180844)['parameters'])
    def parameters(f):
        p=p0.copy();p[8]=(1-f)*p0[8]+f*p1[8];return p
    all_reports=[]
    branches=[(.72,2,None),(.905,4,.9),(.945,8,.94),(.95,16,.945)]
    for frac,order,ancestor in branches:
        if order not in a.orders:continue
        row=rows[frac];p=parameters(frac)
        if order==16:
            info=json.loads(a.long_return.read_text())['accepted_cycles'][0]
        else:info=row['attracting_cycle']
        child_start=np.array(info['refined_start']);child_period=info['period']
        if ancestor is None:
            core_start=np.array(row['primary_refined_start']);core_period=row['primary_period']
        else:
            parent=rows[ancestor]['attracting_cycle']
            tracked=track(np.array(parent['refined_start']),parent['period'],ancestor,frac,parameters,.89)
            if tracked is None:raise RuntimeError(f'Parent shooting failed at {frac}')
            core_start,core_period=tracked
        core_f,core_sol=dynamics(core_start,core_period,p,True)
        child_f,child_sol=dynamics(child_start,child_period,p,True)
        sample_count=65536 if order==16 else 32768
        core=core_sol.sol(np.linspace(0,core_period,sample_count,endpoint=False))[:3].T
        child=child_sol.sol(np.linspace(0,child_period,sample_count,endpoint=False))[:3].T
        scale=np.linalg.norm(np.std(core,axis=0));dist,phase=cKDTree(core).query(np.vstack([child,child[0]]))
        dphase=np.diff(np.unwrap(2*np.pi*phase/len(core)))/(2*np.pi)
        checks=[]
        for n in (640,1280,2560):
            polys=[np.rint(resample(v,n)/scale*1e9).astype(np.int64) for v in (core,child)]
            counts=[]
            for camera in (17,35,81):
                _,mat,_,_=prepare(core[::32],camera)
                val,cross=link_integer(polys[0],polys[1],mat)
                counts.append(dict(matrix=mat.tolist(),link=val,crossings=cross))
            checks.append(dict(vertices=n,projections=counts))
            if n==640:
                views=[dict(camera=c['matrix'],linking=[[0,c['link']],[c['link'],0]],
                            crossings=[[0,c['crossings']],[c['crossings'],0]]) for c in counts]
                (a.out/f'link-{order}.json').write_text(json.dumps(dict(integer_polygons=[v.tolist() for v in polys],views=views),separators=(',',':'))+'\n')
        path=a.out/f'pair-{order}.npz'
        np.savez_compressed(path,points=np.concatenate([core,child]),component=np.repeat([0,1],len(core)),
                            offsets=[0,len(core),2*len(core)])
        parent_poly=dict(family='ring',e=float(p[8]),fraction=frac,child_period_order=order,
                  parameters=p.tolist(),parent_start=core_start.tolist(),child_start=child_start.tolist(),
                  parent=core_f,child=child_f,nearest_phase_turns=float(dphase.sum()),
                  minimum_phase_increment=float(dphase.min()),p95_distance_rms=float(np.quantile(dist,.95)/scale),
                  minimum_distance_rms=float(min(dist)/scale),link_checks=checks,
                  shape=[2*sample_count,3],path=str(path.resolve()),cloud_relative=path.resolve().relative_to(a.out.resolve().parents[1]).as_posix(),
                  point_key='points',component_key='component',point_dtype='float64',sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        all_reports.append(parent_poly)
        (a.out/'report.json').write_text(json.dumps(all_reports,indent=2)+'\n')
        print(order,core_f['dominant'],child_f['dominant'],checks[-1]['projections'][0]['link'],dphase.sum(),flush=True)


if __name__=='__main__':
    main()

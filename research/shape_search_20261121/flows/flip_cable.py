"""Test the observed period-two knot against the primary branch at the same coefficient."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.spatial import cKDTree
from band_continuation import dynamics
from cycle_scan import resample
from exact_knot import exact_diagram,prepare
from symmetry_links import link_integer,check_controls


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',type=Path,default=Path('results/shape_search_20261121/flows'))
    ap.add_argument('--branch',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False);check_controls()
    manifest=json.loads((a.packet/'manifest.json').read_text())
    p0=np.asarray(next(r for r in manifest['stable_cycles'] if r['case']==28756)['parameters'])
    p1=np.asarray(next(r for r in manifest['cases'] if r['case']==180844)['parameters'])
    trials=json.loads((a.branch/'report.json').read_text())['trials']
    row=next(t for t in trials if t.get('cloud_file'));f=row['fraction'];p=p0.copy();p[8]=row['e']
    per=row['primary_period']
    _,sol=dynamics(np.asarray(row['primary_refined_start']),per,p,True)
    core=sol.sol(np.linspace(0,per,32768,endpoint=False))[:3].T
    child=np.load(a.branch/row['cloud_file'])['points']
    scale=np.linalg.norm(np.std(core,axis=0))
    phase_distance,phase=cKDTree(core).query(np.vstack([child,child[0]]),workers=1)
    turns=np.diff(np.unwrap(2*np.pi*phase/len(core)))/(2*np.pi)
    samples=[]
    for size in (640,1280,2560):
        arcs=[resample(c,size) for c in (core,child)]
        ints=[np.rint(c*1e9/scale).astype(np.int64) for c in arcs]
        cams=[];linking=[]
        for k in (17,35,81):
            _,matrix,_,_=prepare(arcs[0],k)
            lk,cross=link_integer(*ints,matrix)
            cams.append(dict(camera=matrix.tolist(),linking=[[0,lk],[lk,0]],crossings=[[0,cross],[cross,0]]))
            linking.append(lk)
        if len(set(linking))>1:raise RuntimeError('Inconsistent linking across cameras')
        samples.append(dict(vertices=size,views=cams))
        if size==640:
            (a.out/'link-witness.json').write_text(json.dumps(dict(integer_polygons=[z.tolist() for z in ints],
                                                  views=cams),separators=(',',':'))+'\n')
    daughter=json.loads((a.branch/row['witness_file']).read_text())
    check=[]
    for size in (640,1280):
        for k in (17,35):
            vertices,mat,_,_=prepare(resample(core,size),k)
            data=exact_diagram(vertices,mat)
            check.append(dict(vertices=size,camera_seed=k,**{key:data[key] for key in ('crossings','coefficients')}))
    tp=np.array(check[0]['coefficients']);cable=np.array(daughter['views'][0]['coefficients'])
    link=abs(samples[-1]['views'][0]['linking'][0][1])
    torus_factor=np.array([(-1)**i for i in range(link)])
    doubled=np.zeros(2*len(tp)-1,dtype=int);doubled[::2]=tp
    product=np.convolve(torus_factor,doubled).tolist()
    report=dict(e=float(p[8]),fraction=f,primary_period=per,child_period=row['attracting_cycle']['period'],
                core_mu=row['primary_multiplier'],child_mu=row['attracting_cycle']['dominant'],scale=scale,
                core_alexander_checks=check,child_alexander=cable.tolist(),linkage=samples,
                projected_phase_turns=float(sum(turns)),worst_backward_jump=float(min(turns)),
                minimum_unsigned_distance_rms=float(min(phase_distance)/scale),
                child_to_core_p95_rms=float(np.quantile(phase_distance,.95)/scale),
                product_polynomial=product,product_matches=product==cable.tolist())
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(a.out/'pair.npz',points=np.concatenate([core,child]),
                        component=np.repeat([0,1],len(core)),offsets=[0,len(core),2*len(core)])
    print(f,link,report['projected_phase_turns'],report['product_matches'],flush=True)


if __name__=='__main__':
    main()

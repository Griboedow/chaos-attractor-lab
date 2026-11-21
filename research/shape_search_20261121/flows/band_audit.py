"""Independent state-only integration and mesh sensitivity for a chaotic band."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
from adaptive import polar_ring


def nearest_segments(points,closed):
    a=closed;b=np.roll(a,-1,axis=0);direction=b-a
    lengths2=(direction*direction).sum(1)
    mid=cKDTree((a+b)/2)
    nearest=cKDTree(a).query(points,workers=1)[0]
    # A minimizing segment is at most the distance to the closest vertex away.
    # Its midpoint therefore lies inside this conservative query ball.
    candidates=mid.query_ball_point(points,nearest+.5*np.sqrt(lengths2.max())+1e-12,workers=1)
    result=np.empty(len(points))
    for i,(x,indices) in enumerate(zip(points,candidates)):
        js=np.asarray(indices)
        v=x-a[js]
        position=np.clip((v*direction[js]).sum(1)/np.maximum(lengths2[js],1e-300),0,1)
        result[i]=np.sqrt(np.min(np.sum((v-position[:,None]*direction[js])**2,axis=1)))
    return result


def summarize(z,orbit,scale):
    zsample=z[::max(1,len(z)//8000)]
    orbit_to=cKDTree(z).query(orbit[::4])[0]/scale
    cloud_to,phase=cKDTree(orbit).query(zsample)
    cloud_to=cloud_to/scale
    counts=np.histogram(phase,np.linspace(0,len(orbit),65))[0]
    q=(.5,.9,.95,.99,1.)
    return dict(size=len(z),rms=float(np.linalg.norm(np.std(z,axis=0))),
                quantiles=list(q),cloud_to_orbit_vertex=np.quantile(cloud_to,q).tolist(),
                orbit_to_cloud=np.quantile(orbit_to,q).tolist(),phase_bins_visited=int(np.sum(counts>0)),
                resolution_segment_p95={str(n):float(np.quantile(nearest_segments(zsample,orbit[::len(orbit)//n]),.95)/scale)
                                        for n in [4096,8192,32768]},
                eight_window_p95=[float(np.quantile(part,.95)) for part in np.array_split(cloud_to,8)])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--case',type=int,default=180844)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    name=f'case-{a.case}.npz'
    with np.load(a.cache/'dense-ring-knot'/name) as f:anchor=f['points'];p=f['parameters']
    orbit=np.load(a.cache/'knot-continuation-a'/name)['points']
    scale=float(np.linalg.norm(np.std(anchor,axis=0)))
    report=dict(case=a.case,common_rms=scale,reference=summarize(anchor,orbit,scale),independent=[])
    for index,method in [(1,'DOP853'),(2,'RK45')]:
        x0=np.load(a.cache/f'controlK{index}'/name)['initial']
        def f(t,y):return polar_ring(y,p)
        warm=solve_ivp(f,(0,900),x0,method=method,rtol=2e-10,atol=2e-12)
        sol=solve_ivp(f,(0,3600),warm.y[:,-1],method=method,rtol=2e-10,atol=2e-12,
                      t_eval=np.linspace(0,3600,72001)[1:])
        if not (warm.success and sol.success):raise RuntimeError('Adaptive integration failed')
        z=sol.y.T
        np.savez_compressed(a.out/f'{method}.npz',points=z)
        report['independent'].append(dict(method=method,rtol=2e-10,atol=2e-12,initial=x0.tolist(),
              burn=900.,observation=3600.,nfev=warm.nfev+sol.nfev,**summarize(z,orbit,scale)))
        (a.out/'geometry.json').write_text(json.dumps(report,indent=2)+'\n')
        print(method,report['independent'][-1]['cloud_to_orbit_vertex'][2],flush=True)


if __name__=='__main__':
    main()

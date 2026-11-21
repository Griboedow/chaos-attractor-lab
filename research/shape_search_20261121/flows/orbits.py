"""Refine short returns and check unstable periodic orbits of the 3D flows."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
from adaptive import full_rhs
from field import velocity_tangent,ring_tangent,rotor_tangent,plane_tangent,crossed_tangent


@njit(cache=True)
def primary_rhs(t,state,p,mode):
    out=np.empty(12)
    column=np.empty(6)
    column[:3]=state[:3]
    for j in range(3):
        column[3:]=state[3+j::3]
        if mode==0:
            values=velocity_tangent(column,p.reshape(3,9))
        elif mode==1:
            values=ring_tangent(column,p)
        elif mode==2:
            values=rotor_tangent(column,p)
        elif mode==3:
            values=plane_tangent(column,p)
        else:
            values=crossed_tangent(column,p)
        out[:3]=values[:3]
        out[3+j::3]=values[3:]
    return out


def seeds(points,dt,limit):
    # Use the least correlated native section with adequate crossings.
    extent=np.sqrt(np.mean(np.sum((points-points.mean(0))**2,axis=1)))
    rows=[]
    for axis in range(3):
        level=np.median(points[:,axis])
        hits=np.flatnonzero((points[:-1,axis]<level)&(points[1:,axis]>=level))
        if len(hits)<3:
            continue
        alpha=(level-points[hits,axis])/(points[hits+1,axis]-points[hits,axis])
        crossings=points[hits]+alpha[:,None]*(points[hits+1]-points[hits])
        times=(hits+alpha)*dt
        for offset in range(1,min(35,len(times))):
            periods=times[offset:]-times[:-offset]
            dist=np.linalg.norm(crossings[offset:]-crossings[:-offset],axis=1)/extent
            good=np.flatnonzero((periods>=1.)&(periods<=limit)&(dist<.06))
            for i in good[np.argsort(dist[good])[:16]]:
                rows.append((float(dist[i])+.0001*periods[i],axis,crossings[i],float(periods[i])))
    return sorted(rows,key=lambda r:r[0]),extent


def refine(axis,start,period,scale,p,mode):
    free=[i for i in range(3) if i!=axis]
    current=start.copy()
    initial_period=period
    for iteration in range(13):
        y=np.r_[current,np.eye(3).ravel()]
        sol=solve_ivp(primary_rhs,(0.,period),y,args=(p,mode),method='DOP853',
                      rtol=3e-11,atol=3e-13,t_eval=[period])
        if not sol.success:
            return None
        end=sol.y[:,-1]
        delta=end[:3]-current
        if np.linalg.norm(delta)/scale<1.5e-9:
            return current,period
        matrix=end[3:].reshape(3,3)-np.eye(3)
        velocity=primary_rhs(0.,end,p,mode)[:3]
        system=np.column_stack([matrix[:,free],velocity])
        try:
            change=np.linalg.solve(system,-delta)
        except np.linalg.LinAlgError:
            return None
        damping=max(1.,np.linalg.norm(change[:2])/(.18*scale),abs(change[2])/(.1*initial_period))
        current[free]+=change[:2]/damping
        period+=change[2]/damping
        if not .7*initial_period < period < 1.3*initial_period:
            return None
    return None


def check(start,period,scale,tree,p,mode):
    grid=np.linspace(0,period,401)
    y=np.r_[start,np.eye(3).ravel(),0.]
    sol=solve_ivp(full_rhs,(0.,period),y,args=(p,mode),method='DOP853',rtol=4e-12,atol=4e-14,
                  t_eval=grid,dense_output=True)
    if not sol.success:
        return None
    cloud=sol.y[:3].T
    err=np.linalg.norm(cloud[-1]-start)/scale
    if min(np.linalg.norm(sol.sol(period/k)[:3]-start)/scale for k in range(2,13))<1e-5:
        return None
    eig=np.linalg.eigvals(sol.y[3:12,-1].reshape(3,3))
    one=np.argmin(abs(eig-1))
    rest=np.delete(eig,one)
    large=rest[np.argmax(abs(rest))]
    if err>1e-6 or abs(eig[one]-1)>1e-3 or abs(large)<1.02 or abs(large.imag)>1e-4:
        return None
    dist=tree.query(cloud,workers=1)[0]/scale
    if np.quantile(dist,.95)>.08 or np.linalg.norm(np.ptp(cloud,axis=0))/scale<.12:
        return None
    return dict(initial_xyz=start.tolist(),period=float(period),independent_closure_rms_units=float(err),
                neutral_multiplier=float(eig[one].real),unstable_multiplier=float(large.real),
                stable_multiplier_from_liouville=float(np.exp(sol.y[12,-1])/large.real),
                integrated_divergence=float(sol.y[12,-1]),cloud_distance_p50_rms_units=float(np.median(dist)),
                cloud_distance_p95_rms_units=float(np.quantile(dist,.95)),
                orbit_span_rms_units=float(np.linalg.norm(np.ptp(cloud,axis=0))/scale))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True,type=Path)
    ap.add_argument('--family',choices=('cartesian','ring','rotor','planes','crossed'),required=True)
    ap.add_argument('--ids',nargs='+',type=int,required=True)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--limit',type=float,default=85.)
    ap.add_argument('--attempts',type=int,default=30)
    ap.add_argument('--retain',type=int,default=3)
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    mode=('cartesian','ring','rotor','planes','crossed').index(args.family)
    for case in args.ids:
        with np.load(args.source/f'case-{case}.npz') as f:
            points=f['points']
            p=f['parameters'].ravel()
            dt=float(f['step_size']*f['sample_stride'])
        candidates,scale=seeds(points,dt,args.limit)
        tree=cKDTree(points)
        results=[]; attempted=[]; guesses=[]
        for score,axis,x,period in candidates:
            if any(abs(period-value)<.015*period for value in attempted):
                continue
            if any(a==axis and abs(period-t)<.03*period and np.linalg.norm(x-y)<.12*scale for a,y,t in guesses):
                continue
            guesses.append((axis,x,period))
            got=refine(axis,x,period,scale,p,mode)
            if got is not None:
                report=check(*got,scale,tree,p,mode)
                if report and not any(abs(got[1]-v)<1e-4*got[1] for v in attempted):
                    results.append(report)
                    attempted.append(got[1])
            if len(results)>=args.retain or len(guesses)>=args.attempts:
                break
        (args.out/f'case-{case}.json').write_text(json.dumps(dict(case=case,family=args.family,
            guesses=len(guesses),orbits=results),indent=2)+'\n')
        print(case,len(results),flush=True)


if __name__=='__main__':
    main()

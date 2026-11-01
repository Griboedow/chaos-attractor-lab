"""Search common full-state return times after minimizing circular symmetry."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
import spatial_fast as m
import symmetry


def series(t,q,k):
    q=q.copy()
    for j in range(1,len(q)):
        q[j,np.sum(q[j-1]*q[j],axis=1)<0]*=-1
    qfun=CubicSpline(t,q);kfun=CubicSpline(t,k)
    def state(time):
        qq=qfun(time);qq/=np.linalg.norm(qq,axis=-1,keepdims=True)
        return qq,kfun(time)
    return state


def search(t,q,k,scale,lower=.2,upper=5.0,samples=241):
    assert upper>lower>0 and t[-1]-t[0]>upper+5
    state=series(t,q,k)
    anchors=np.linspace(max(t[0]+upper,t[-1]-12),t[-1],7)
    reference=[state(v) for v in anchors]
    def objective(lag):
        values=[]
        for time,(qa,ka) in zip(anchors,reference):
            qb,kb=state(time-lag)
            values.append(symmetry.distance(qa,ka,qb,kb,scale)['rms']**2)
        return float(np.mean(values))
    lags=np.linspace(lower,upper,samples)
    scores=np.array([objective(lag) for lag in lags])
    indices=np.flatnonzero((scores[1:-1]<scores[:-2])&(scores[1:-1]<scores[2:]))+1
    minima=[]
    for j in indices:
        fit=minimize_scalar(objective,bounds=(lags[j-1],lags[j+1]),method='bounded',options={'xatol':1e-10})
        lag=float(fit.x);reduced=[];ordinary=[]
        for time,(qa,ka) in zip(anchors,reference):
            qb,kb=state(time-lag)
            reduced.append(symmetry.distance(qa,ka,qb,kb,scale))
            ordinary.append(symmetry.distance(qa,ka,qb,kb,scale,reduce=False))
        minima.append(dict(lag=lag,rms_over_anchors=float(np.sqrt(fit.fun)),
            reduced_rms_max=max(v['rms'] for v in reduced),
            reduced_orientation_max=max(v['orientation_max'] for v in reduced),
            reduced_curvature_max=max(v['curvature_max'] for v in reduced),
            alignment_angles=[v['angle'] for v in reduced],
            ordinary_rms_max=max(v['rms'] for v in ordinary),optimizer_success=bool(fit.success)))
    return dict(anchors=anchors.tolist(),search_lags=lags.tolist(),search_squared_scores=scores.tolist(),
        minima=minima,best=min(minima,key=lambda x:x['rms_over_anchors']) if minima else None,
        scope='Finite lag scan and scalar numerical minimization using seven full-state pairs. A low residual supports a numerical relative recurrence; it does not prove a period, irrational rotation, attraction, chaos or continuum convergence.')


def check():
    t=np.linspace(0,20,2001);n=5;period=1.7;omega=.731
    q=[];k=[]
    for time in t:
        phase=2*np.pi*time/period
        qq=np.array([m.expq(np.array([.2*(j+1)/n*np.cos(phase),.13*np.sin(phase),.03*np.cos(phase)])) for j in range(n)])
        kk=np.array([[j+np.cos(phase),.4*np.sin(phase)] for j in range(n)])
        qq,kk=symmetry.action(qq,kk,omega*time);q.append(qq);k.append(kk)
    result=search(t,np.array(q),np.array(k),5,lower=1.3,upper=2.1)
    best=result['best']
    assert abs(best['lag']-period)<1e-6 and best['reduced_rms_max']<1e-7
    assert best['ordinary_rms_max']>.1
    return dict(analytic_shape_period=period,analytic_angular_drift=omega*period,recovered=best)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('folders',nargs='*',type=Path)
    ap.add_argument('--check',action='store_true')
    ap.add_argument('--upper',type=float,default=5.0)
    ap.add_argument('--samples',type=int,default=241)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    if args.check:
        result=check()
    else:
        result=[]
        for folder in args.folders:
            report=json.loads((folder/'report.json').read_text(encoding='utf-8'))
            assert report['parameters']['eta']==1
            with np.load(folder/'trajectory.npz') as d:
                if d['t'][-1]-d['t'][0]<12:continue
                row=search(d['t'],d['q'],d['k'],abs(report['parameters']['chi']),upper=args.upper,samples=args.samples)
            row['run']=folder.name;row['success']=report['success'];result.append(row)
            print(json.dumps(dict(run=folder.name,best=row['best'])),flush=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if args.check:print(json.dumps(result))

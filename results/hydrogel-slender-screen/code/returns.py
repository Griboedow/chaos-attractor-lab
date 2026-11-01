"""Compare circular-rod full states on a rotation-invariant tip-height section."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import symmetry


def analyze(folder):
    report=json.loads((folder/'report.json').read_text(encoding='utf-8'))
    assert report['parameters']['eta']==1
    with np.load(folder/'trajectory.npz') as data:
        t=data['t'];q=data['q'];k=data['k'];tip=data['tip']
    # Make the interpolated quaternion representative continuous in time.
    for j in range(1,len(q)):
        signs=np.sum(q[j-1]*q[j],axis=1)<0
        q[j,signs]*=-1
    cutoff=max(t[-1]/2,t[-1]-20)
    mask=t>=cutoff
    level=float(np.mean(tip[mask,2]))
    z=CubicSpline(t,tip[:,2]-level)
    qfun=CubicSpline(t,q);kfun=CubicSpline(t,k)
    ids=np.flatnonzero((tip[:-1,2]<level)&(tip[1:,2]>=level)&(t[:-1]>=cutoff))
    roots=np.array([brentq(z,t[j],t[j+1],xtol=1e-13) for j in ids])
    qq=qfun(roots);qq/=np.linalg.norm(qq,axis=-1,keepdims=True);kk=kfun(roots)
    scale=abs(report['parameters']['chi'])
    returns=[]
    for lag in range(1,min(9,len(roots))):
        ordinary=[];reduced=[]
        for j in range(max(lag,len(roots)-10),len(roots)):
            args=(qq[j-lag],kk[j-lag],qq[j],kk[j],scale)
            ordinary.append(symmetry.distance(*args,reduce=False))
            reduced.append(symmetry.distance(*args,reduce=True))
        returns.append(dict(lag=lag,pairs=len(reduced),
            ordinary_rms_max=max(x['rms'] for x in ordinary),
            reduced_rms_max=max(x['rms'] for x in reduced),
            reduced_rms_last=reduced[-1]['rms'],
            reduced_rms_values=[x['rms'] for x in reduced],
            reduced_orientation_max=max(x['orientation_max'] for x in reduced),
            reduced_curvature_max=max(x['curvature_max'] for x in reduced),
            alignment_angles=[x['angle'] for x in reduced]))
    best=min(returns,key=lambda r:r['reduced_rms_max']) if returns else None
    probes=np.linspace(cutoff,t[-1],21)
    pq=qfun(probes);pq/=np.linalg.norm(pq,axis=-1,keepdims=True);pk=kfun(probes)
    spread=[symmetry.distance(q[-1],k[-1],a,b,scale)['rms'] for a,b in zip(pq,pk)]
    result=dict(run=folder.name,success=report['success'],completed_time=float(t[-1]),
        statistics_start=float(cutoff),section='Upward tip-height crossing of the late-window mean',
        section_height=level,section_times=roots.tolist(),periods=np.diff(roots).tolist(),
        height_range=np.ptp(tip[mask,2]).item(),
        tip_radius_range=np.ptp(np.linalg.norm(tip[mask,:2],axis=1)).item(),
        curvature_scale=scale,return_errors=returns,best_tested_lag=best,
        quotient_spread_from_final=dict(times=probes.tolist(),rms=spread,maximum=max(spread)),
        scope='Finite late-window returns of cubic-interpolated full orientation/curvature states; SO(2) alignment is exact for this distance, but no asymptotic, chaos or continuum conclusion follows.')
    print(json.dumps(dict(run=folder.name,success=report['success'],events=len(roots),best=best)))
    return result


if __name__=='__main__':
    root=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('folders',nargs='*',type=Path)
    ap.add_argument('--output',type=Path,default=root/'returns.json')
    args=ap.parse_args()
    folders=args.folders or sorted(root.iterdir())
    rows=[analyze(p) for p in folders if (p/'report.json').exists() and (p/'trajectory.npz').exists()]
    args.output.write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')

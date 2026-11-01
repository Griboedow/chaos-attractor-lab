import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import spatial_fast as m

root=Path(__file__).resolve().parent
summaries=[]
for folder in sorted(root.iterdir()):
    if not folder.is_dir() or not (folder/'report.json').exists() or not (folder/'trajectory.npz').exists(): continue
    r=json.loads((folder/'report.json').read_text(encoding='utf-8'))
    d=np.load(folder/'trajectory.npz');t=d['t'];tip=d['tip'];q=d['q'];k=d['k']
    cutoff=max(t[-1]/2,t[-1]-20)
    mask=t>=cutoff
    axis=1 if r['plane']=='yz' else 0
    spline=CubicSpline(t,tip[:,axis]);qfun=CubicSpline(t,q);kfun=CubicSpline(t,k)
    ids=np.flatnonzero((tip[:-1,axis]<0)&(tip[1:,axis]>=0)&(t[:-1]>=cutoff))
    roots=np.array([brentq(spline,t[j],t[j+1],xtol=1e-13) for j in ids])
    returns=[]
    if len(roots)>1:
        qq=qfun(roots);qq/=np.linalg.norm(qq,axis=-1,keepdims=True);kk=kfun(roots)
        for lag in range(1,min(9,len(roots))):
            angles=[]
            for a,b in zip(qq[:-lag],qq[lag:]):
                angles.append(max(np.linalg.norm(m.logq(m.mul(np.r_[x[0],-x[1:]],y))) for x,y in zip(a,b)))
            returns.append(dict(lag=lag,orientation_max=float(max(angles[-10:])),
                                curvature_max=float(np.max(abs(kk[lag:]-kk[:-lag])[-10:]))))
    xy=tip[mask,:2]-tip[mask,:2].mean(axis=0)
    sv=np.linalg.svd(xy,compute_uv=False)
    finalgrad=m.energy_gradient(q[-1],k[-1],m.parameters(**dict(lam=r['parameters']['slenderness'],eta=r['parameters']['eta'],xi=r['parameters']['xi'],chi=r['parameters']['chi'],nu=r['parameters']['nu']))[0])[1]
    out=dict(run=folder.name,success=r['success'],duration=float(t[-1]),statistics_start=float(cutoff),
             xy_minor_major=float(sv[1]/sv[0]) if sv[0]>0 else 0,
             periods=np.diff(roots).tolist(),return_errors=returns,
             mean_period=float(np.mean(np.diff(roots))) if len(roots)>1 else None,
             period_std=float(np.std(np.diff(roots))) if len(roots)>1 else None,
             final_tip_sample_speed=float(np.linalg.norm(tip[-1]-tip[-2])/(t[-1]-t[-2])),
             final_elastic_gradient_max=float(np.max(abs(finalgrad))))
    summaries.append(out)
    print(json.dumps({key:out[key] for key in ['run','success','mean_period','period_std','xy_minor_major','final_tip_sample_speed']}))
(root/'analysis.json').write_text(json.dumps(summaries,indent=2)+'\n',encoding='utf-8')

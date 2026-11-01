import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
root=Path(__file__).resolve().parent
ref=np.load('D:/Codex/hydrogel-filament-20261101/event-16-BDF/trajectory.npz')
ref_fun=CubicSpline(ref['t'],ref['y'].T)
rows=[]
for name in ['planar-16-dt002','planar-16-dt001']:
    d=np.load(root/name/'trajectory.npz');t=d['t'];q=d['q'];k=d['k']
    state=np.c_[-2*np.arctan2(q[:,:,1],q[:,:,0]),-k[:,:,0]]
    expected=ref_fun(t)
    spline=CubicSpline(t,d['tip'][:,1])
    ids=np.flatnonzero((d['tip'][:-1,1]<0)&(d['tip'][1:,1]>=0)&(t[:-1]>15))
    roots=[brentq(spline,t[i],t[i+1],xtol=1e-13) for i in ids]
    rows.append(dict(run=name,max_angle_difference=float(np.max(abs(state[:,:16]-expected[:,:16]))),
                     max_curvature_difference=float(np.max(abs(state[:,16:]-expected[:,16:]))),
                     mean_period_after15=float(np.mean(np.diff(roots))),
                     periods_after15=np.diff(roots).tolist()))
result=dict(reference='Earlier adaptive BDF planar angle solver: 16 segments, rtol=1e-7, atol=1e-9; cubic interpolation of saved dt=0.02 states.',rows=rows,
            angle_error_ratio=rows[0]['max_angle_difference']/rows[1]['max_angle_difference'],
            curvature_error_ratio=rows[0]['max_curvature_difference']/rows[1]['max_curvature_difference'])
(root/'planar-comparison.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))

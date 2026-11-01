import json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parent
rows=[]
for name in ['rotation-n16','rotation-n32','rotation-n64']:
    d=np.load(root/name/'trajectory.npz');mask=d['t']>=d['t'][-1]-20
    t=d['t'][mask];tip=d['tip'][mask];radius=np.linalg.norm(tip[:,:2],axis=1)
    angle=np.unwrap(np.arctan2(tip[:,1],tip[:,0]))
    fit=np.polyfit(t-t.mean(),angle,1)
    row=dict(run=name,time_window=[float(t[0]),float(t[-1])],angular_rate=float(fit[0]),
             inferred_period=float(2*np.pi/abs(fit[0])),angle_fit_rms=float(np.sqrt(np.mean((angle-np.polyval(fit,t-t.mean()))**2))),
             radius_mean=float(radius.mean()),radius_min=float(radius.min()),radius_max=float(radius.max()),
             height_min=float(tip[:,2].min()),height_max=float(tip[:,2].max()))
    rows.append(row);print(json.dumps(row))
(root/'rotation-summary.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')

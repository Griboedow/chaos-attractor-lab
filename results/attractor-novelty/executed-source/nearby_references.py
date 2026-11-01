from pathlib import Path
import json
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
base=Path(__file__).parent/'reference'

@njit
def henon(m1,m2,b):
 z=np.array([.1,.1,.1]);points=np.zeros((16000,3))
 for i in range(56000):
  x,y,w=z;z[0]=y;z[1]=w;z[2]=m1+b*x+m2*y-w*w
  if i>=40000:points[i-40000]=z
 return points

records=[]
for m1,m2,b in [(0.,.85,.7),(0.,.815,.7),(1.4,.2,.1),(-.02,.85,.7)]:
 name=f'Henon3D-{m1}-{m2}-{b}'
 points=henon(m1,m2,b);assert np.isfinite(points).all()
 np.savez_compressed(base/'trajectories'/f'{name}.npz',points=points)
 records.append({'name':name,'map':'(y,z,M1+B*x+M2*y-z*z)','M1':m1,'M2':m2,'B':b,'discarded':40000,'observed':16000,'start':[.1,.1,.1],
 'sources':['https://www.ma.imperial.ac.uk/~dturaev/mypapers/bc2005.pdf','https://doi.org/10.1063/5.0123426']})
for a in [1.1,1.2,1.3,1.5,1.6,1.7,1.8,1.85,1.89,1.95,2.0]:
 def rhs(t,s):
  x,y,z=s
  return [-a*x-4*y-4*z-y*y,-a*y-4*z-4*x-z*z,-a*z-4*x-4*y-x*x]
 name=f'Halvorsen-a{a}'
 sol=solve_ivp(rhs,(0,1200),[-6.4,0,0],method='DOP853',rtol=1e-9,atol=1e-11,t_eval=np.linspace(400,1200,16000))
 row={'name':name,'a':a,'b':4,'span':[0,1200],'sample_span':[400,1200],'start':[-6.4,0,0],'success':bool(sol.success),'n':16000,'method':'DOP853','rtol':1e-9,'atol':1e-11}
 if sol.success and sol.y.shape==(3,16000):np.savez_compressed(base/'trajectories'/f'{name}.npz',points=sol.y.T)
 records.append(row);print(name,sol.success,flush=True)
(base/'nearby-runs.json').write_text(json.dumps(records,indent=2)+'\n')

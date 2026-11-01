"""Independent transcription of the four external atlas equations."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--out',type=Path,required=True)
base=parser.parse_args().out
(base/'trajectories').mkdir(parents=True,exist_ok=True)

def rhs(t,s,kind):
 x,y,z=s;r=x*x+y*y
 if kind=='Aurelia':
  return [1.4*(z-.5)*x-1.9*y+.5*(x*x-y*y),1.9*x+1.4*(z-.5)*y-x*y,1.4*(1+z)-z**3-1.9*r]
 if kind=='Naiad':
  return [(1.2*(z-.68)+1.7)*x-3.7*y,3.7*x+(1.2*(z-.68)-1.7)*y,1.2+2*z-z**3-2.1*r]
 if kind=='Cassiopea':
  return [2*(z-.85)*x-1.8*y+.9*(x**3-3*x*y*y),1.8*x+2*(z-.85)*y+.9*(y**3-3*x*x*y),1.28+1.8*z-z**3-2.8*r]
 return [(1.4104*z*z-.6701+.5534*z)*x-2.4266*y+.3*(x**3-3*x*y*y),
         (1.4104*z*z-.6701-.5534*z)*y+2.4266*x+.3*(y**3-3*x*x*y),
         3.3*z-z**3-2.1355*z*r+2*1.2153*x*y]

rows=[]
for kind in ['Aurelia','Naiad','Cassiopea','Mobula']:
 sol=solve_ivp(lambda t,z:rhs(t,z,kind),(0,2000),[.3,.2,.1],method='DOP853',rtol=1e-10,atol=1e-12,t_eval=np.linspace(400,2000,16000))
 assert sol.success and sol.y.shape==(3,16000)
 np.savez_compressed(base/'trajectories'/f'Atlas-{kind}.npz',points=sol.y.T)
 rows.append({'name':kind,'source':'https://github.com/Sakeeb91/attractor-atlas/blob/821ebb816d3e90bb9f589556f48cb455c4d9f545/README.md','span':[0,2000],'sample_span':[400,2000],'n':16000,'rtol':1e-10,'atol':1e-12,'success':bool(sol.success)})
 print(kind,flush=True)
(base/'atlas-runs.json').write_text(json.dumps(rows,indent=2)+'\n')

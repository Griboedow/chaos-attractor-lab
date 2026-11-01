"""Independent low-order coefficient checks and endpoint convergence comparisons."""
from pathlib import Path
import contextlib,io,json,runpy,sys
import numpy as np
from flint import arb,ctx

base=Path(__file__).parent;out=base/'icosa-taylor-checks';out.mkdir(exist_ok=False)
sys.argv=['icosa_taylor_precision.py','--benchmark']
with contextlib.redirect_stdout(io.StringIO()):ns=runpy.run_path(str(base/'icosa_taylor_precision.py'))
advance=ns['advance'];exact=ns['exact'];d=np.load(base/'icosa-v1/input.npz');c=d['centers'][177];u=d['axes'][177];s=d['sigma'][177];g=d['gamma'][177]
skew=np.zeros((60,3,3));skew[:,0,1]=-u[:,2];skew[:,0,2]=u[:,1];skew[:,1,0]=u[:,2];skew[:,1,2]=-u[:,0];skew[:,2,0]=-u[:,1];skew[:,2,1]=u[:,0]
checks=[]
for x in [np.array(ns['row']['point']),np.array([.31,-.82,.4]),np.array([1.,.2,-.3])]:
    r=x-c;v=np.cross(u,r);w=np.exp(-.5*(r*r).sum(1)/s**2);f=-g*x+(w[:,None]*v).sum(0);J=-g*np.eye(3)+(w[:,None,None]*(skew-v[:,:,None]*r[:,None,:]/s**2)).sum(0)
    a=np.array([float(v) for v in advance([exact(v) for v in x],arb(1),1)[0]])
    b=np.array([float(v) for v in advance([exact(v) for v in x],arb(1),2)[0]])
    checks.append({'point':x.tolist(),'first_order_error':float(np.linalg.norm(a-(x+f))),'second_order_error':float(np.linalg.norm(b-(x+f+.5*J@f)))})
assert max(max(r['first_order_error'],r['second_order_error']) for r in checks)<1e-12
ctx.dps=120;rows=[]
def read(index,digits,order,h):return json.loads((base/f'icosa-taylor-precision/{index}-full-d{digits}-o{order}-h{h}.json').read_text())
def distance(a,b):return sum((arb(x)-arb(y))**2 for x,y in zip(a['endpoint'],b['endpoint'])).sqrt().str(12,radius=False)
for index in [0,20]:
    a,b,c0,e=read(index,80,24,.04),read(index,80,24,.02),read(index,80,32,.04),read(index,110,32,.04)
    rows.append({'index':index,'order24_step_halving_difference':distance(a,b),'order24_h002_vs_order32_h004':distance(b,c0),'order32_80_vs_110_digits':distance(c0,e),'high_precision_closure_error':e['closure_error']})
report={'coefficient_checks':checks,'endpoint_checks':rows,'limits':'Consistency and convergence evidence only. Truncation error is not enclosed and the binary64 shooting initial data are not exact periodic roots.'}
(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)

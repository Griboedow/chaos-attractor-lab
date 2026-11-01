"""Independent implementation, event derivative and step refinement controls."""
from pathlib import Path
from return_paths import INPUT,STRUCTURE,CYCLES,OUTPUT
import json,time
import numpy as np
import icosa_poincare as flow
from icosa_fast_return import rhs,rk_return

base=OUTPUT;out=base/'icosa-return-controls';out.mkdir(exist_ok=False)
rng=np.random.default_rng(2026110157);error=0.
for _ in range(30):
    x=np.r_[rng.normal(size=3),rng.normal(size=6)];error=max(error,float(np.max(abs(rhs(0,x)-flow.field(0,x)))))
assert error<1e-11
flow.field=rhs
old=json.loads((base/'icosa-return-strips/report.json').read_text());selected=[]
for strip in [0,1]:
    rows=[r for r in old['rows'] if r['strip']==strip]
    selected.extend([rows[i] for i in [0,2,13,24,26]])
checks=[];begin=time.time()
for row in selected:
    uv=np.array(row['input']);q,t,J,speed=flow.evaluate(uv,rtol=3e-14,atol=3e-16,max_step=.015);fixed=[]
    for h in [.005,.0025,.00125]:
        a,b,c=rk_return(uv,h);fixed.append({'h':h,'output':a.tolist(),'time':float(b),'output_difference':float(np.linalg.norm(a-q)),'time_difference':float(abs(b-t))})
    checks.append({'strip':row['strip'],'input':uv.tolist(),'output':q.tolist(),'time':t,'derivative':J.tolist(),'speed':speed,'adaptive_output_difference':float(np.linalg.norm(q-row['output'])),'adaptive_derivative_difference':float(np.max(abs(J-np.array(row['derivative'])))),'fixed_steps':fixed})
    print(json.dumps({'strip':row['strip'],'input':uv.tolist(),'fine_rk4_error':fixed[-1]['output_difference'],'seconds':time.time()-begin}),flush=True)
# A central difference tests the event-time correction, including the stable input.
derivatives=[]
for uv in [np.array([0.,0.]),np.array([.0021125,0.])]:
    q,t,J,speed=flow.evaluate(uv,rtol=3e-14,atol=3e-16,max_step=.015);eps=1e-8
    approx=np.column_stack([(flow.evaluate(uv+eps*e,rtol=3e-14,atol=3e-16,max_step=.015)[0]-flow.evaluate(uv-eps*e,rtol=3e-14,atol=3e-16,max_step=.015)[0])/(2*eps) for e in np.eye(2)])
    derivatives.append({'input':uv.tolist(),'analytic':J.tolist(),'central_difference':approx.tolist(),'max_absolute_error':float(np.max(abs(J-approx)))})
report={'implementation_max_error':error,'rows':checks,'derivative_controls':derivatives,'seconds':time.time()-begin,'limits':'Independent numerical checks only. No interval enclosure of a whole strip or exact first-return certification.'};(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps({'implementation_max_error':error,'derivative_errors':[v['max_absolute_error'] for v in derivatives]}),flush=True)

from paths import OUTPUT, ROOT, RESEARCH
"""Longer validation of every positive from the extended survey."""
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import brentq
from anderson_extended import survey,step,g_and_derivative

root=OUTPUT
out=root/'anderson-long';out.mkdir(exist_ok=False)
d=np.load(root/'anderson-v2/data.npz');ids=d['positive_ids']
# Independent coefficient-space solve and complex-step differentiation.
def direct(z,p):
    k,c,reg,kind=p
    if kind==4:g=c+k*z
    elif kind==2:g=c+k*np.tanh(z)
    elif kind==3:g=c+k*z/(1+z*z)
    elif kind==0:g=c+k*np.sin(z)
    else:g=c+k*np.cos(z)
    r=g-z
    coef=np.linalg.solve(np.outer(r,r)+reg*np.eye(3),np.ones(3))
    coef/=coef.sum()
    return np.array([z[1],z[2],np.dot(coef,g)])

rng=np.random.default_rng(2026110134);eqerr=0.;jacerr=0.
for i in range(100):
    z=rng.normal(size=3);p=np.array([rng.uniform(.5,8),rng.uniform(-2,2),10**rng.uniform(-2,2),[0,2,3,4][i%4]])
    actual,_=step(z,np.ones(3),p);expected=direct(z,p)
    eqerr=max(eqerr,float(np.max(np.abs(actual-expected))))
    for j in range(3):
        v=np.eye(3)[j];_,deriv=step(z,v,p)
        independent=np.imag(direct(z.astype(complex)+1e-25j*v,p))/1e-25
        jacerr=max(jacerr,float(np.max(np.abs(deriv-independent))))
assert eqerr<1e-8 and jacerr<1e-7,(eqerr,jacerr)
# Every contraction pilot has a unique root; independently bracket it.
old=np.load(root/'anderson-v1/data.npz');rooterror=0.;residual=0.
for p,z in zip(old['parameters'],old['end']):
    def residual_function(x):return g_and_derivative(x,p[0],p[1],int(p[3]))[0]-x
    fixed=brentq(residual_function,-20,20,xtol=1e-14)
    rooterror=max(rooterror,float(np.max(np.abs(z-fixed))));residual=max(residual,abs(residual_function(z[-1])))
checks={'direct_KKT_max_error':eqerr,'complex_step_jacobian_max_error':jacerr,'contraction_pilot_max_root_error':rooterror,'contraction_pilot_max_residual':residual,'contraction_pilot_status_counts':{str(int(v)):int((old['status']==v).sum()) for v in np.unique(old['status'])}}
(out/'checks.json').write_text(json.dumps(checks,indent=2))
begin=time.time()
w,tail,end,status,steps,bounds=survey(d['parameters'][ids],d['starts'][ids],50000,500000,20000)
np.savez_compressed(out/'data.npz',ids=ids,parameters=d['parameters'][ids],starts=d['starts'][ids],windows=w,tail=tail,end=end,status=status,steps=steps,max_abs=bounds)
print(json.dumps({'seconds':time.time()-begin,'checks':checks,'positive_ids':ids[(status==1)&(w.min(1)>.005)].tolist(),'counts':{str(int(v)):int((status==v).sum()) for v in np.unique(status)}}),flush=True)

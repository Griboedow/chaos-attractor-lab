"""Independent step/Jacobian checks and checks of apparent convergence."""
from pathlib import Path
import json,hashlib
import numpy as np
from follow import map_jac,integrate
from pilot import run
root=Path('D:/Codex/unexpected-attractors-20261101/fastica')
def reference(w,X,kind,alpha):
 u=X@w
 if kind==0:g=np.tanh(alpha*u);gp=alpha*(1-g*g)
 elif kind==1:e=np.exp(-u*u/2);g=u*e;gp=(1-u*u)*e
 else:g=u**3;gp=3*u*u
 h=X.T@g/len(X)-gp.mean()*w
 return h/np.linalg.norm(h)
rng=np.random.default_rng(2026110103);err=0.;jacerr=0.;epsilon=2e-6;rows=[]
for case in [15,40,204,246,351]:
 d=json.loads((root/f'follow-01/case-{case}-input.json').read_text());X=np.array(d['X']);kind=d['kind'];alpha=d['alpha']
 for i in range(64):
  w=rng.normal(size=4);w/=np.linalg.norm(w);v=rng.normal(size=4);v-=w*(v@w);v/=np.linalg.norm(v);f,J,hn=map_jac(w,X,kind,alpha)
  err=max(err,float(abs(f-reference(w,X,kind,alpha)).max()))
  a=w+epsilon*v;b=w-epsilon*v;a/=np.linalg.norm(a);b/=np.linalg.norm(b);fd=(reference(a,X,kind,alpha)-reference(b,X,kind,alpha))/(2*epsilon)
  jacerr=max(jacerr,float(np.linalg.norm(fd-J@v)/max(1,np.linalg.norm(J@v))))
 if case==40:
  w=np.load(root/'long-01/case-40.npz')['endpoint'];f,J,hn=map_jac(w,X,kind,alpha);signed=np.linalg.norm(f-w);anti=np.linalg.norm(f+w);sign=1 if signed<anti else -1;Q=np.linalg.svd(w.reshape(1,4),full_matrices=True)[2][1:].T
  rows.append({'case':40,'endpoint':w.tolist(),'one_step_residual':float(signed),'antipodal_residual':float(anti),'return_multiplier_moduli':sorted(abs(np.linalg.eigvals(Q.T@(sign*J)@Q)).tolist()),'meaning':'Numerical projective fixed point; this is not an interval certificate.'})
 if case==246:
  lib=json.loads((root/'library-controls-v2.json').read_text());w=np.array(next(r['terminal_weight'] for r in lib if r['case']==246));f,J,hn=map_jac(w,X,kind,alpha);result=run(X,w,kind,alpha,1000,1000000,100000)
  rows.append({'case':246,'library_stopping_state':w.tolist(),'one_step_projective_distance':float(min(np.linalg.norm(f-w),np.linalg.norm(f+w))),'continuation_warmup':1000,'continuation_observed':1000000,'continuation_rates':result[0].tolist(),'period_up_to_sign':int(result[2]),'status':int(result[4]),'meaning':'Continue the recorded library stopping state without the tolerance-based stopping rule.'})
assert err<1e-12 and jacerr<2e-7
report={'seed':2026110103,'states_per_case':64,'maximum_step_error':err,'maximum_tangent_relative_error':jacerr,'finite_difference_step':epsilon,'outcomes':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(root/'controls.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

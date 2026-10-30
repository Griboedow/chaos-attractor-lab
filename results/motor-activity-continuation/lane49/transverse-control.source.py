import sys,json,hashlib,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment
sys.path.insert(0,str(Path('research').resolve()))
import motor_postcritical_cycles as cyc
import motor_postcritical_branch as branch
p=Path('results/motor-activity-continuation/lane49/periodic/cycle.npz');a=np.load(p);x=a['quotient_initial_state'];pa=a['parameters'];T=float(a['time'][-1]);w=np.ones(4)/4
s=cyc.integrate(x,T,pa,w,dense=True,rtol=1e-12,max_step=.0125)
blocks=[]
for j in range(4):
 def fun(t,b):
  c,ss,xx=cyc.unpack(s.sol(t),w);v=-pa[0]*xx+pa[1]*(ss-ss.mean());r=pa[1];rho=pa[0]
  mat=np.array([[-1,-v[j]-r*ss[j],rho*ss[j]],[v[j],-1+r*c[j],-rho*c[j]],[0,r,-rho]])
  return (mat@b.reshape(3,3)).ravel()
 answers={}
 for method in ['DOP853','Radau']:
  z=solve_ivp(fun,(0,T),np.eye(3).ravel(),method=method,rtol=1e-12,atol=1e-14,max_step=.01);assert z.success;answers[method]=z.y[:,-1].reshape(3,3)
 eig=np.linalg.eigvals(answers['DOP853']);blocks.append(dict(cluster=j,real=eig.real.tolist(),imag=eig.imag.tolist(),modulus=abs(eig).tolist(),DOP853_Radau_matrix_error=float(np.max(abs(answers['DOP853']-answers['Radau'])))))
def rows(q):return np.stack(cyc.unpack(q,w),axis=1)
y=rows(x);z=rows(s.sol(T/3));cost=np.linalg.norm(z[:,None]-y[None,:],axis=2);ri,ci=linear_sum_assignment(cost)
record=dict(status='complete',source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),period=T,transverse_blocks=blocks,third_period_permutation=ci.tolist(),third_period_state_error=float(cost[ri,ci].max()),scope='Numerical independent contrast-block solver checks; dense mean trajectory is DOP853, not an interval enclosure.')
out=Path('results/motor-activity-continuation/lane49/transverse-control.json');assert not out.exists();out.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');print(record)

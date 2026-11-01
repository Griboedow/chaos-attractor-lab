import json
import time
from pathlib import Path
import numpy as np
import spatial as slow
import spatial_fast as fast
root=Path(__file__).resolve().parent
rng=np.random.default_rng(260409487);n=8
beta,zeta,xi,chi,eta=fast.parameters()
q=np.array([fast.expq(v) for v in rng.normal(scale=0.1,size=(n,3))])
qold=fast.rotate_increment(q,rng.normal(scale=0.003,size=(n,3)))
qprev=fast.rotate_increment(qold,rng.normal(scale=0.003,size=(n,3)))
no=fast.geometry(qold)[1];npv=fast.geometry(qprev)[1]
ko=rng.normal(scale=0.5,size=(n,2));kp=ko+rng.normal(scale=0.01,size=(n,2));dt=0.002
checks=[]
for order in (1,2):
    jac=fast.mechanical_jacobian(q,ko,kp,no,npv,beta,zeta,xi,chi,eta,dt,order)
    fd=np.empty_like(jac)
    for j in range(3*n):
        inc=np.zeros((n,3));inc[j//3,j%3]=1e-6
        plus=fast.residual_q(fast.rotate_increment(q,inc),ko,no,kp,npv,beta,zeta,xi,chi,eta,dt,order)[0]
        minus=fast.residual_q(fast.rotate_increment(q,-inc),ko,no,kp,npv,beta,zeta,xi,chi,eta,dt,order)[0]
        fd[:,j]=((plus-minus)/(2e-6)).ravel()
    error=float(np.max(abs(fd-jac))/np.max(abs(fd)))
    assert error<1e-8,error
    checks.append(dict(order=order,jacobian_relative_error=error))
q,k,_=fast.consistent_initial(8,beta,'spatial');nodes=fast.geometry(q)[1]
args=(q,k,nodes,k.copy(),nodes.copy(),beta,zeta,xi,chi,eta,dt,1,np.zeros((8,3)),1e-8)
a=slow.implicit_step(*args);b=fast.implicit_step(*args)
assert a[4] and b[4]
error=max(float(np.max(abs(a[0]-b[0]))),float(np.max(abs(a[1]-b[1]))))
assert error<1e-10,error
timing={}
for name,module in [('full_difference',slow),('local_jacobian',fast)]:
    start=time.perf_counter()
    for _ in range(30): module.implicit_step(*args)
    timing[name]=(time.perf_counter()-start)/30
result=dict(jacobian_checks=checks,one_step_max_difference=error,seconds_per_step=timing)
(root/'fast-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))

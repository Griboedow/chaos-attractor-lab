"""Independent finite differences and earlier-assembly controls for faster drag."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import spatial_fast as old
import spatial_accelerated as new

root=Path(__file__).resolve().parent
rng=np.random.default_rng(260409489)
rows=[]
for eta in (.73,1.0):
    n=9
    q=np.array([old.expq(v) for v in rng.normal(scale=.2,size=(n,3))])
    k=rng.normal(scale=.8,size=(n,2));kp=k+rng.normal(scale=.03,size=k.shape)
    no=old.geometry(old.rotate_increment(q,rng.normal(scale=.002,size=(n,3))))[1]
    npv=old.geometry(old.rotate_increment(q,rng.normal(scale=.002,size=(n,3))))[1]
    beta,zeta,xi,chi,eta=old.parameters(lam=127,eta=eta,chi=117)
    for order in (1,2):
        a=old.drag_jacobian(q,no,npv,zeta,xi,.006,order)
        b=new.drag_jacobian(q,no,npv,zeta,xi,.006,order)
        fd=np.empty_like(b)
        for j in range(3*n):
            inc=np.zeros((n,3));inc[j//3,j%3]=1e-6
            plus=old.drag_gradient(old.rotate_increment(q,inc),no,npv,zeta,xi,.006,order)
            minus=old.drag_gradient(old.rotate_increment(q,-inc),no,npv,zeta,xi,.006,order)
            fd[:,j]=((plus-minus)/2e-6).ravel()
        assembly=float(np.max(abs(a-b))/np.max(abs(a)))
        independent=float(np.max(abs(fd-b))/np.max(abs(fd)))
        full_a=old.mechanical_jacobian(q,k,kp,no,npv,beta,zeta,xi,chi,eta,.006,order)
        full_b=new.mechanical_jacobian(q,k,kp,no,npv,beta,zeta,xi,chi,eta,.006,order)
        full=float(np.max(abs(full_a-full_b))/np.max(abs(full_a)))
        assert assembly<1e-12 and independent<1e-8 and full<1e-12
        rows.append(dict(eta=eta,order=order,drag_assembly_relative_error=assembly,
            drag_finite_difference_relative_error=independent,full_jacobian_relative_error=full))
seed=Path('D:/Codex/hydrogel-slender-search-20261101/lambda127-field117/trajectory.npz')
with np.load(seed) as d:q=d['q'][-1];k=d['k'][-1]
beta,zeta,xi,chi,eta=old.parameters(lam=127,eta=1,chi=117)
no=old.geometry(q)[1]
args=(q,k,no,k,no,beta,zeta,xi,chi,eta,.006,1,np.zeros((len(q),3)),1e-8)
a=old.implicit_step(*args);b=new.implicit_step(*args)
assert a[4] and b[4]
step_error=max(float(np.max(abs(a[i]-b[i]))) for i in range(3))
assert step_error<1e-10
times={}
for name,module in [('previous',old),('accelerated',new)]:
    start=time.perf_counter()
    for repeat in range(30):module.mechanical_jacobian(q,k,k,no,no,beta,zeta,xi,chi,eta,.006,2)
    times[name]=(time.perf_counter()-start)/30
result=dict(checks=rows,one_step_maximum_difference=step_error,jacobian_seconds=times,
    timing_scope='Warm 64-segment Jacobian calls on this machine while other runs are active; not a universal speed claim.',
    source_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ['spatial_fast.py','spatial_accelerated.py','check_accelerated.py']})
(root/'accelerated-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))

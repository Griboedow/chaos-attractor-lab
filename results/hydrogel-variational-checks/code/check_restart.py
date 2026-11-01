"""Check interpolation locations, twist projection and saved BDF history."""
import hashlib
import json
from pathlib import Path
import numpy as np
import restart as driver
import spatial_fast as m
import symmetry

root=Path(__file__).resolve().parent
seed=Path('D:/Codex/hydrogel-slender-search-20261101/lambda127-field117/trajectory.npz')
with np.load(seed) as d:q=d['q'][-1];k=d['k'][-1]
beta,zeta,xi,chi,eta=m.parameters(lam=127,eta=1,chi=117)
same_q,same_k=driver.remesh(q,k,64)
assert np.array_equal(same_q,q) and np.array_equal(same_k,k)
qr,kr=driver.remesh(q,k,128)
projected,projection=driver.consistent_twist(qr,kr,beta)
centerline_error=float(np.max(abs(m.geometry(qr)[1]-m.geometry(projected)[1])))
assert centerline_error<1e-13 and projection['final_scaled_torque']<1e-10
angle=.641
qa,ka=symmetry.action(q,k,angle)
ar,ak=driver.remesh(qa,ka,128)
br,bk=symmetry.action(qr,kr,angle)
remesh_equivariance=max(float(np.max(abs(ar-br))),float(np.max(abs(ak-bk))))
assert remesh_equivariance<1e-12
# A constant-strain director field is reproduced by the interpolation and
# its short endpoint extrapolation. This control is independent of the seed.
n=11;v=np.array([.2,-.1,.07]);s=(np.arange(n)+.5)/n
qtest=np.array([m.expq(t*v) for t in s]);ktest=np.zeros((n,2))
actual,_=driver.remesh(qtest,ktest,23)
exact=np.array([m.expq(t*v) for t in (np.arange(23)+.5)/23])
constant_strain_error=float(np.max(abs(actual-exact)))
assert constant_strain_error<1e-14
# Replay the pilot's first 20 accepted steps, then compare a 21st step from
# live arrays against the same step from its written checkpoint.
state_q=projected.copy();state_k=kr.copy();nodes=m.geometry(state_q)[1]
kp=state_k.copy();npv=nodes.copy();predictor=np.zeros((128,3))
for step in range(1,21):
    qn,kn,inc,its,ok,res,hinge=m.implicit_step(state_q,state_k,nodes,kp,npv,beta,zeta,xi,chi,eta,.006,1 if step==1 else 2,predictor,1e-8)
    assert ok
    kp,npv=state_k.copy(),nodes.copy();state_q,state_k=qn,kn
    nodes=m.geometry(state_q)[1];predictor=inc
with np.load(root/'pilot128/checkpoint.npz') as c:
    fields=['q','k','nodes','k_previous','nodes_previous','predictor']
    live=[state_q,state_k,nodes,kp,npv,predictor]
    saved=[c[name].copy() for name in fields]
    checkpoint_errors={name:float(np.max(abs(a-b))) for name,a,b in zip(fields,live,saved)}
    assert all(e==0 for e in checkpoint_errors.values())
    a=m.implicit_step(state_q,state_k,nodes,kp,npv,beta,zeta,xi,chi,eta,.006,2,predictor,1e-8)
    b=m.implicit_step(saved[0],saved[1],saved[2],saved[3],saved[4],beta,zeta,xi,chi,eta,.006,2,saved[5],1e-8)
    assert a[4] and b[4]
    next_error=max(float(np.max(abs(a[i]-b[i]))) for i in range(3))
    assert next_error==0
result=dict(identity_remesh_exact=True,constant_strain_interpolation_maximum_error=constant_strain_error,
    remesh_equivariance_maximum_error=remesh_equivariance,twist_projection=projection,
    twist_projection_centerline_maximum_change=centerline_error,
    replayed_pilot_checkpoint_errors=checkpoint_errors,loaded_checkpoint_next_step_error=next_error,
    seed_sha256=hashlib.sha256(seed.read_bytes()).hexdigest(),
    scope='Initial-state interpolation, algebraic consistency and exact saved multistep-history checks; no asymptotic or spatial-convergence certificate.')
(root/'restart-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))

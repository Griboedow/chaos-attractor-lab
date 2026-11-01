"""Check a complete tangent step at the 128-segment field-124 state."""
import json
from pathlib import Path
import numpy as np
import spatial_accelerated as m
import tangent

root=Path(__file__).resolve().parent
with np.load(root/'field124-n128-relax/checkpoint.npz') as c:
    q=c['q'];k=c['k'];nodes=c['nodes'];kp=c['k_previous'];npv=c['nodes_previous'];predictor=c['predictor']
n=len(q);beta,zeta,xi,chi,eta=m.parameters(lam=127,eta=1,chi=124)
args=(beta,zeta,xi,chi,eta,.006,2,predictor,1e-11)
result=m.implicit_step(q,k,nodes,kp,npv,*args);assert result[4]
qn,kn=result[:2];ez=np.array([0.,0.,1.])
gk=np.stack([-k[:,1],k[:,0]],axis=1)[:,:,None]
gkp=np.stack([-kp[:,1],kp[:,0]],axis=1)[:,:,None]
gn=np.cross(ez,nodes)[:,:,None];gnp=np.cross(ez,npv)[:,:,None]
answer=tangent.propagate(qn,k,kp,nodes,npv,gk,gkp,gn,gnp,beta,zeta,xi,chi,eta,.006,2)
expected=[np.array([m.matrix(a).T@ez-ez for a in qn]),
    np.stack([-kn[:,1],kn[:,0]],axis=1),np.cross(ez,m.geometry(qn)[1])]
symmetry_errors=[float(np.max(abs(a[:,:,0]-b))) for a,b in zip(answer,expected)]
assert max(symmetry_errors)<1e-6,symmetry_errors
rng=np.random.default_rng(260409493)
theta=rng.normal(size=(n,3,1));dk=rng.normal(size=(n,2,1))*chi
dn=tangent.node_variation(q,theta);zerosk=np.zeros_like(dk);zerosn=np.zeros_like(dn)
answer=tangent.propagate(qn,k,kp,nodes,npv,dk,zerosk,dn,zerosn,beta,zeta,xi,chi,eta,.006,2)
eps=1e-6;outputs=[]
for sign in (1,-1):
    qt=m.rotate_increment(q,sign*eps*theta[:,:,0])
    value=m.implicit_step(qt,k+sign*eps*dk[:,:,0],m.geometry(qt)[1],kp,npv,*args)
    assert value[4];outputs.append(value)
angles=[np.array([m.logq(m.mul(np.r_[a[0],-a[1:]],b)) for a,b in zip(qn,v[0])]) for v in outputs]
observed=[(angles[0]-angles[1])/(2*eps),(outputs[0][1]-outputs[1][1])/(2*eps),
    (m.geometry(outputs[0][0])[1]-m.geometry(outputs[1][0])[1])/(2*eps)]
errors=[float(np.max(abs(a[:,:,0]-b))/max(1.,np.max(abs(b)))) for a,b in zip(answer,observed)]
assert max(errors)<1e-5,errors
out=dict(segments=n,chi=chi,neutral_symmetry_errors=symmetry_errors,
    complete_step_difference_errors=errors,finite_difference_step=eps,
    scope='One refined-state step, not a long spectrum or convergence proof.')
(root/'refined-tangent-check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out))

"""Validate BDF tangent propagation against complete-step differences and symmetry."""
import json
from pathlib import Path
import numpy as np
import spatial_accelerated as m
import tangent

root=Path(__file__).resolve().parent
rng=np.random.default_rng(260409491);n=12;rows=[]
q=np.array([m.expq(v) for v in rng.normal(scale=.04,size=(n,3))])
qprev=m.rotate_increment(q,rng.normal(scale=.002,size=(n,3)))
k=rng.normal(scale=.5,size=(n,2));kp=k+rng.normal(scale=.02,size=(n,2))
nodes=m.geometry(q)[1];npv=m.geometry(qprev)[1]
theta=rng.normal(size=(n,3,1));thetap=rng.normal(size=(n,3,1))
dk=rng.normal(size=(n,2,1));dkp=rng.normal(size=(n,2,1))
dn=tangent.node_variation(q,theta);dnp=tangent.node_variation(qprev,thetap)
beta,zeta,xi,chi,eta=m.parameters(lam=127,eta=1,chi=117)
for order in (1,2):
    args=(beta,zeta,xi,chi,eta,.006,order,np.zeros((n,3)),1e-11)
    result=m.implicit_step(q,k,nodes,kp,npv,*args);assert result[4]
    qn,kn=result[:2]
    exact=tangent.propagate(qn,k,kp,nodes,npv,dk,dkp,dn,dnp,beta,zeta,xi,chi,eta,.006,order)
    eps=1e-5
    outputs=[]
    for sign in (1,-1):
        qt=m.rotate_increment(q,sign*eps*theta[:,:,0]);qp=m.rotate_increment(qprev,sign*eps*thetap[:,:,0])
        value=m.implicit_step(qt,k+sign*eps*dk[:,:,0],m.geometry(qt)[1],kp+sign*eps*dkp[:,:,0],m.geometry(qp)[1],*args)
        assert value[4];outputs.append(value)
    angles=[]
    for value in outputs:
        angles.append(np.array([m.logq(m.mul(np.r_[a[0],-a[1:]],b)) for a,b in zip(qn,value[0])]))
    observed=[(angles[0]-angles[1])/(2*eps),(outputs[0][1]-outputs[1][1])/(2*eps),
        (m.geometry(outputs[0][0])[1]-m.geometry(outputs[1][0])[1])/(2*eps)]
    errors=[float(np.max(abs(a[:,:,0]-b))/max(1.,np.max(abs(b)))) for a,b in zip(exact,observed)]
    assert max(errors)<1e-6,errors
    ez=np.array([0.,0.,1.])
    gk=np.stack([-k[:,1],k[:,0]],axis=1)[:,:,None]
    gkp=np.stack([-kp[:,1],kp[:,0]],axis=1)[:,:,None]
    gn=np.cross(ez,nodes)[:,:,None];gnp=np.cross(ez,npv)[:,:,None]
    predicted=tangent.propagate(qn,k,kp,nodes,npv,gk,gkp,gn,gnp,beta,zeta,xi,chi,eta,.006,order)
    expected_theta=np.array([m.matrix(a).T@ez-ez for a in qn])
    expected_k=np.stack([-kn[:,1],kn[:,0]],axis=1)
    expected_n=np.cross(ez,m.geometry(qn)[1])
    symmetry_errors=[float(np.max(abs(a[:,:,0]-b))) for a,b in zip(predicted,[expected_theta,expected_k,expected_n])]
    assert max(symmetry_errors)<1e-7,symmetry_errors
    rows.append(dict(order=order,complete_step_difference_errors=errors,neutral_symmetry_errors=symmetry_errors))
out=dict(checks=rows,finite_difference_step=1e-5,
    scope='Two BDF orders checked against independently perturbed nonlinear solves and the exact circular rotation tangent. Numerical checks, not interval bounds.')
(root/'tangent-check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out))

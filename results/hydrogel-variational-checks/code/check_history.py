"""Check that QR rescaling transforms both time levels of the BDF tangent."""
import json
from pathlib import Path
import numpy as np
import spatial_accelerated as m
import tangent

root=Path(__file__).resolve().parent
rng=np.random.default_rng(260409492);n=12;p=3
beta,zeta,xi,chi,eta=m.parameters(lam=127,eta=1,chi=117)
q,k,_=m.consistent_initial(n,beta,'spatial');nodes=m.geometry(q)[1]
theta=rng.normal(scale=.1,size=(n,3,p));dk=rng.normal(size=(n,2,p))
dn=tangent.node_variation(q,theta)
args=(beta,zeta,xi,chi,eta,.006)
first=m.implicit_step(q,k,nodes,k,nodes,*args,1,np.zeros((n,3)),1e-11)
assert first[4];q1,k1=first[:2];n1=m.geometry(q1)[1]
a=tangent.propagate(q1,k,k,nodes,nodes,dk,dk,dn,dn,*args,1)
vector=np.concatenate([a[0].reshape(3*n,p),a[1].reshape(2*n,p)/abs(chi)],axis=0)/np.sqrt(n)
basis,upper=np.linalg.qr(vector,mode='reduced');inverse=np.linalg.solve(upper,np.eye(p))
transform=lambda v:np.einsum('nip,pq->niq',v,inverse)
second=m.implicit_step(q1,k1,n1,k,nodes,*args,2,first[2],1e-11)
assert second[4]
direct=tangent.propagate(second[0],k1,k,n1,nodes,a[1],dk,a[2],dn,*args,2)
scaled=tangent.propagate(second[0],k1,k,n1,nodes,transform(a[1]),transform(dk),transform(a[2]),transform(dn),*args,2)
errors=[float(np.max(abs(np.einsum('nip,pq->niq',b,upper)-a))/max(1.,np.max(abs(a)))) for a,b in zip(direct,scaled)]
assert max(errors)<1e-8,errors
out=dict(qr_orthogonality_error=float(np.max(abs(basis.T@basis-np.eye(p)))),
    two_step_basis_change_relative_errors=errors,
    scope='A full three-direction basis change is propagated through both current and previous BDF states. This verifies the normalization algebra, not asymptotic exponent convergence.')
(root/'history-check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out))

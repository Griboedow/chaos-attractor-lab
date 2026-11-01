"""Check the conditional-probability factorization of the instantaneous Jacobian."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.special import logsumexp
import model,model_dual

def verify(L,z,wx,wy):
    m,n=L.shape;p=m-1;d=m+n-2;la,lb=model.targets(L)
    zz=model_dual.step(z,np.ones(d),L,la,lb,wx,wy)[0]
    J=np.column_stack([model_dual.step(z,v,L,la,lb,wx,wy)[1] for v in np.eye(d)])
    joint=L+np.r_[0.,zz[:p]][:,None]+np.r_[0.,z[p:]][None,:]
    logs=joint-(logsumexp(joint,axis=1)[:,None]+logsumexp(joint,axis=0)[None,:])/2
    s=np.linalg.svd(np.exp(logs),compute_uv=False);mu=s[1:]**2
    rx,ry=1-wx,1-wy;expected=[]
    for b in mu:expected.extend(np.roots([1,-(rx+ry+wx*wy*b),rx*ry]).astype(complex))
    expected.extend([complex(rx)]*max(0,m-n));expected.extend([complex(ry)]*max(0,n-m))
    actual=np.linalg.eigvals(J).astype(complex);cost=abs(actual[:,None]-np.array(expected)[None,:]);a,b=linear_sum_assignment(cost)
    error=float(cost[a,b].max());assert error<1e-6
    determinant=float(abs(np.linalg.det(J)-rx**(m-1)*ry**(n-1)));assert determinant<1e-10
    state_error=None
    if wx==wy:
        state_error=float(np.max(abs(zz-model.step(z,np.ones(d),L,la,lb,wx)[0])));assert state_error==0
    return dict(shape=[m,n],wx=wx,wy=wy,singular_values=s.tolist(),spectral_radius=float(max(abs(actual))),
        expected_eigenvalues=[[v.real,v.imag] for v in expected],maximum_eigenvalue_matching_error=error,determinant_error=determinant,equal_rate_state_error=state_error)

rows=[];rng=np.random.default_rng(201212810)
for m,n in [(2,2),(2,3),(3,2),(2,4),(4,2),(3,3),(4,4),(3,5)]:
    for wx,wy in [(.8,.8),(1.5,1.5),(1.99,1.99),(.3,1.95),(1.95,.3),(2.3,.7)]:
        for _ in range(5):
            L=rng.normal(size=(m,n))*2;z=rng.normal(size=m+n-2)*3;r=verify(L,z,wx,wy)
            if max(wx,wy)<2:assert r['spectral_radius']<1+1e-10
            rows.append(r)
result=dict(checks=rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Numerical consistency checks of an algebraic factorization; sampled eigenvalues are not a proof of the universal spectral condition.')
Path('schur-check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(cases=len(rows),max_eigen_error=max(r['maximum_eigenvalue_matching_error'] for r in rows),max_det_error=max(r['determinant_error'] for r in rows))))

"""Independent probability-domain update, finite differences and determinant checks."""
from pathlib import Path
import json
import numpy as np
from scipy.special import logsumexp
import model

rng=np.random.default_rng(201212563);rows=[]
for m,n in [(2,2),(3,4),(4,3)]:
    L=rng.normal(size=(m,n));la,lb=model.targets(L);d=m+n-2
    K=np.exp(L);a=K.sum(axis=1);b=K.sum(axis=0)
    for omega in [.8,1.,1.7,1.99,2.,2.1]:
        z=rng.normal(scale=.3,size=d);v=rng.normal(size=d)
        actual,directional=model.step(z,v,L,la,lb,omega)
        u=np.exp(np.r_[0,z[:m-1]]);vv=np.exp(np.r_[0,z[m-1:]])
        un=u**(1-omega)*(a/(K@vv))**omega
        vn=vv**(1-omega)*(b/(K.T@un))**omega
        independent=np.r_[np.log(un[1:]/un[0]),np.log(vn[1:]/vn[0])]
        update=float(np.max(abs(actual-independent)))
        eps=1e-6
        fd=(model.step(z+eps*v,v,L,la,lb,omega)[0]-model.step(z-eps*v,v,L,la,lb,omega)[0])/(2*eps)
        derivative=float(np.max(abs(fd-directional)))
        jac=model.jacobian(z,L,la,lb,omega)
        determinant=float(abs(np.linalg.det(jac)-(1-omega)**d))
        assert update<1e-12 and derivative<1e-7 and determinant<1e-11
        fixed=float(np.max(abs(model.step(np.zeros(d),v,L,la,lb,omega)[0])))
        assert fixed<1e-12
        row=dict(shape=[m,n],omega=omega,update_error=update,directional_error=derivative,determinant_error=determinant,fixed_point_error=fixed)
        if omega<2:
            bound=model.bounds(L,la,lb,omega)
            corners=rng.uniform(-1,1,size=(100,d))*bound
            excess=max(float(np.max(abs(model.step(p,v,L,la,lb,omega)[0])-bound)) for p in corners)
            row['sampled_invariant_box_excess']=excess
            assert excess<1e-12
        rows.append(row)
result=dict(checks=rows,scope='Floating-point implementation checks; the determinant and box formulas are derived algebraically, not certified by the random samples.')
(Path(__file__).resolve().parent/'model-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(cases=len(rows),maximum_update_error=max(r['update_error'] for r in rows),maximum_directional_error=max(r['directional_error'] for r in rows),maximum_determinant_error=max(r['determinant_error'] for r in rows))))

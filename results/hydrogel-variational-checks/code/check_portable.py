"""Synthetic assembly and remeshing controls without external trajectory inputs."""
import json
from pathlib import Path
import numpy as np
import spatial_fast as old
import spatial_accelerated as new
import restart
import symmetry

rng=np.random.default_rng(260409493);rows=[]
for eta in (.73,1.):
    q=np.array([old.expq(v) for v in rng.normal(scale=.1,size=(9,3))])
    no=old.geometry(old.rotate_increment(q,rng.normal(scale=.001,size=(9,3))))[1]
    npv=old.geometry(old.rotate_increment(q,rng.normal(scale=.001,size=(9,3))))[1]
    beta,zeta,xi,chi,eta=old.parameters(lam=127,eta=eta,chi=124)
    for order in (1,2):
        a=old.drag_jacobian(q,no,npv,zeta,xi,.006,order)
        b=new.drag_jacobian(q,no,npv,zeta,xi,.006,order)
        direction=rng.normal(size=(9,3));eps=1e-6
        plus=old.drag_gradient(old.rotate_increment(q,eps*direction),no,npv,zeta,xi,.006,order)
        minus=old.drag_gradient(old.rotate_increment(q,-eps*direction),no,npv,zeta,xi,.006,order)
        fd=((plus-minus)/(2*eps)).ravel()
        error=float(np.max(abs(b@direction.ravel()-fd))/np.max(abs(fd)))
        assembly=float(np.max(abs(a-b))/np.max(abs(a)))
        assert assembly<1e-12 and error<1e-7
        rows.append(dict(eta=eta,order=order,assembly_relative_error=assembly,directional_difference_relative_error=error))
n=11;v=np.array([.2,-.1,.07]);q=np.array([old.expq(s*v) for s in (np.arange(n)+.5)/n]);k=np.zeros((n,2))
qr,kr=restart.remesh(q,k,23)
expected=np.array([old.expq(s*v) for s in (np.arange(23)+.5)/23])
interpolation=float(np.max(abs(qr-expected)));assert interpolation<1e-14
qa,ka=symmetry.action(q,k,.641);ar,ak=restart.remesh(qa,ka,23);br,bk=symmetry.action(qr,kr,.641)
equivariance=max(float(np.max(abs(ar-br))),float(np.max(abs(ak-bk))));assert equivariance<1e-12
out=dict(assembly=rows,constant_strain_interpolation_error=interpolation,remesh_equivariance_error=equivariance)
(Path(__file__).resolve().parent/'portable-check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out))

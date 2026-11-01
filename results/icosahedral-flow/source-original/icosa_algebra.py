"""Check numerical representations of identities for the ideal equivariant field."""
from pathlib import Path
import json
import numpy as np
from tetra_spectra import field

root=Path(__file__).parent;d=np.load(root/'icosa-v1/input.npz');out=root/'icosa-algebra';out.mkdir(exist_ok=False)
group=d['group'];closure=max(min(np.linalg.norm(a@b-c) for c in group) for a in group for b in group)
orthogonality=float(np.max(np.abs(group@group.transpose(0,2,1)-np.eye(3))));deterror=float(np.max(np.abs(np.linalg.det(group)-1)))
rows=[];rng=np.random.default_rng(2026110151)
for ident in [177,183,43,29,85]:
    c=d['centers'][ident];u=d['axes'][ident];s=d['sigma'][ident];g=d['gamma'][ident]
    def fn(z):
        r=z-c;return -g*z+(np.exp(-.5*(r*r).sum(1)/s**2)[:,None]*np.cross(u,r)).sum(0)
    f0,j0=field(np.zeros(3),c,u,s,g);ferr=jerr=diverr=equiv=0.
    for z in rng.normal(size=(12,3)):
        f,j=field(z,c,u,s,g);jc=np.column_stack([np.imag(fn(z+1e-25j*e))/1e-25 for e in np.eye(3)])
        ferr=max(ferr,float(np.max(np.abs(f-fn(z)))));jerr=max(jerr,float(np.max(np.abs(j-jc))));diverr=max(diverr,abs(float(np.trace(j)+3*g)))
        for a in group:equiv=max(equiv,float(np.max(np.abs(fn(a@z)-a@fn(z)))))
    assert max(ferr,jerr,diverr,equiv)<1e-10
    row={'id':ident,'seed_center':d['seed_center'][ident].tolist(),'seed_axis':d['seed_axis'][ident].tolist(),'sigma':float(s),'gamma':float(g),'amplitude':float(d['amplitude'][ident]),'start':d['starts'][ident].tolist(),'origin_field_residual':float(np.linalg.norm(f0)),'origin_jacobian_scalar_residual':float(np.max(np.abs(j0+g*np.eye(3)))),'field_implementation_error':ferr,'complex_step_jacobian_error':jerr,'divergence_error':diverr,'equivariance_error':equiv,'absorbing_radius_threshold':float(np.linalg.norm(u,axis=1).sum()*s*np.exp(-.5)/g)}
    rows.append(row)
(out/'report.json').write_text(json.dumps({'group_size':len(group),'closure_error':closure,'orthogonality_error':orthogonality,'determinant_error':deterror,'cases':rows,'limits':'Exact symmetry and the exact origin statement refer to the ideal algebraic rotation group. Stored binary64 rotations approximate that group; their residuals are reported here. Divergence and the radial bound hold for arbitrary stored centers and axes as real coefficients.'},indent=2))
print(json.dumps(rows),flush=True)

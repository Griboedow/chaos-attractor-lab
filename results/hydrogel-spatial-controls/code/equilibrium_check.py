"""Construct a spurious planar fixed point of the coarse midpoint scheme."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import spatial_fast as m

n=16;alpha=33.0;h=1/n
theta=np.zeros(n)
theta[1]=brentq(lambda x:x-alpha*h*np.sin(x/2),1e-7,np.pi,xtol=1e-14)
for i in range(2,n):
    prev=theta[i-1]
    theta[i]=brentq(lambda x:x-prev-alpha*h*np.sin((prev+x)/2),prev-alpha*h,prev+alpha*h,xtol=1e-14)
k=np.zeros((n,2))
k[1:,0]=-np.diff(theta)/h
q=np.array([m.expq(np.array([-v,0.0,0.0])) for v in theta])
beta,zeta,xi,chi,eta=m.parameters(eta=1,chi=-alpha)
energy,gradient,hinge=m.energy_gradient(q,k,beta)
drive=m.target_curvature(q,chi,eta)
nodes=m.geometry(q)[1]
residual=m.residual_q(q,k,nodes,k,nodes,beta,zeta,xi,chi,eta,0.004,2)[0]
assert np.max(abs(residual))<1e-10
assert np.max(abs(drive-k))<1e-10
result=dict(segments=n,chi=-alpha,eta=1,theta=theta.tolist(),curvature=k.tolist(),
            tip=nodes[-1].tolist(),energy=float(energy),scaled_mechanical_residual=float(np.max(abs(residual))),
            relaxation_residual=float(np.max(abs(drive-k))),max_hinge_angle=float(hinge),
            interior_feedback_ratio=alpha/(2*n),first_cell_feedback_ratio=alpha/(4*n),
            scope='Constructed discrete fixed point, not an equilibrium of the continuum source model. Nonzero first interior root exists when alpha/(2N)>1; the clamp half-cell is still unique when alpha/(4N)<1.')
observed=np.load(Path(__file__).resolve().parent/'rotation-n16/trajectory.npz')
tip=observed['tip'][-1]
azimuth=np.arctan2(-tip[0],tip[1]);z=m.expq(np.array([0.0,0.0,azimuth]));zc=z.copy();zc[1:]*=-1
target=np.array([m.mul(m.mul(z,v),zc) for v in q])
qobs=observed['q'][-1]
errors=[np.linalg.norm(m.logq(m.mul(np.r_[a[0],-a[1:]],b))) for a,b in zip(target,qobs)]
ktarget=(m.matrix(z)@np.c_[k,np.zeros(n)].T).T[:,:2]
result['comparison_at_time80']=dict(azimuth=float(azimuth),max_orientation_angle=float(max(errors)),
                                  max_curvature_difference=float(np.max(abs(ktarget-observed['k'][-1]))),
                                  tip_distance=float(np.linalg.norm(m.matrix(z)@nodes[-1]-tip)))
(Path(__file__).resolve().parent/'spurious-equilibrium.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({key:result[key] for key in ['tip','energy','scaled_mechanical_residual','relaxation_residual','interior_feedback_ratio']}))

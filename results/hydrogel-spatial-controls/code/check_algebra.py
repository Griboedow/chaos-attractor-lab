import json
import sys
from pathlib import Path
import numpy as np
import spatial as m
sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')
import hydrogel_planar as planar

root=Path(__file__).resolve().parent
rng=np.random.default_rng(260409486)
n=10
beta,zeta,xi,chi,eta=m.parameters()
q=np.array([m.expq(v) for v in rng.normal(scale=0.12,size=(n,3))])
k=rng.normal(scale=0.3,size=(n,2))
energy,gradient,_=m.energy_gradient(q,k,beta)
fd=np.empty_like(gradient)
for i in range(n):
    for j in range(3):
        inc=np.zeros((n,3));inc[i,j]=1e-6
        ep=m.energy_gradient(m.rotate_increment(q,inc),k,beta)[0]
        em=m.energy_gradient(m.rotate_increment(q,-inc),k,beta)[0]
        fd[i,j]=(ep-em)/2e-6
energy_error=np.max(abs(fd-gradient))/np.max(abs(fd))
assert energy_error<1e-8,energy_error

qold=m.rotate_increment(q,rng.normal(scale=1e-3,size=(n,3)))
frames,nodes=m.geometry(q)
oldnodes=m.geometry(qold)[1]
dt=0.003
analytic=m.drag_gradient(q,oldnodes,oldnodes,zeta,xi,dt,1)
quad=np.zeros_like(analytic)
points,weights=np.polynomial.legendre.leggauss(3)
for segment in range(n):
    drag=(frames[segment]*zeta)@frames[segment].T*xi
    for f,w in zip((points+1)/2,weights/2):
        velocity=((1-f)*(nodes[segment]-oldnodes[segment])+f*(nodes[segment+1]-oldnodes[segment+1]))/dt
        force=drag@velocity
        for joint in range(segment+1):
            local=frames[joint].T@force
            quad[joint]+=w/(n*n)*(f if joint==segment else 1)*np.array([-local[1],local[0],0])
drag_error=np.max(abs(quad-analytic))/np.max(abs(quad))
assert drag_error<1e-12,drag_error

theta=rng.normal(scale=0.1,size=n)
curvature=rng.normal(scale=0.3,size=n)
velocity=rng.normal(scale=0.1,size=n)
qplanar=np.array([m.expq(np.array([-v,0.0,0.0])) for v in theta])
kplanar=np.c_[-curvature,np.zeros(n)]
grad=m.energy_gradient(qplanar,kplanar,beta)[1]
p=planar.parameters()
prhs=planar.angle_rhs_raw(np.r_[theta,curvature],*planar.rhs_parameters(p))
mass=planar.angle_mass(theta,p['normal_drag'],p['tangent_drag'],xi)
planar_energy_error=np.max(abs(mass@prhs[:n]-grad[:,0]))/np.max(abs(grad[:,0]))
assert planar_energy_error<1e-10,planar_energy_error
drive=m.target_curvature(qplanar,chi,eta)
planar_curvature_error=np.max(abs(-drive[:,0]-curvature-prhs[n:]))
assert planar_curvature_error<1e-12,planar_curvature_error
perturbed=np.array([m.expq(np.array([-v,0.0,0.0])) for v in theta+1e-7*velocity])
g=m.drag_gradient(perturbed,m.geometry(qplanar)[1],m.geometry(qplanar)[1],zeta,xi,1e-7,1)
planar_drag_error=np.max(abs(-g[:,0]-mass@velocity))/np.max(abs(mass@velocity))
assert planar_drag_error<1e-7,planar_drag_error
assert np.max(abs(grad[:,1:]))<1e-20
assert np.max(abs(drive[:,1]))<1e-20

result=dict(energy_gradient_relative_error=float(energy_error),drag_quadrature_relative_error=float(drag_error),
            planar_energy_relative_error=float(planar_energy_error),planar_curvature_absolute_error=float(planar_curvature_error),
            planar_drag_relative_error=float(planar_drag_error),random_seed=260409486,segments=n)
(root/'algebra.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))

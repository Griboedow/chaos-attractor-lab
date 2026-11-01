"""SO(2) comparison for a circular-section rod, using complete material states."""
from pathlib import Path
import json
import numpy as np
import spatial_fast as model


def action(q,k,angle):
    """R -> Z R Z^T, kappa -> Z kappa; the clamped identity stays fixed."""
    c,s=np.cos(angle),np.sin(angle)
    qq=q.copy();kk=k.copy()
    qq[...,1]=c*q[...,1]-s*q[...,2]
    qq[...,2]=s*q[...,1]+c*q[...,2]
    kk[...,0]=c*k[...,0]-s*k[...,1]
    kk[...,1]=s*k[...,0]+c*k[...,1]
    return qq,kk


def frames(q):
    return np.array([model.matrix(v) for v in q])


def best_angle(qa,ka,qb,kb,curvature_scale):
    """Global minimum of frame Frobenius + scaled curvature distance.

    Conjugating a frame gives only first and second angular harmonics.
    Every stationary angle is obtained from the tan(angle/2) quartic.
    Quaternion signs therefore have no effect on the comparison.
    """
    a,b=frames(qa),frames(qb)
    va=a[:,0,2]+1j*a[:,1,2];vb=b[:,0,2]+1j*b[:,1,2]
    wa=a[:,2,0]+1j*a[:,2,1];wb=b[:,2,0]+1j*b[:,2,1]
    ua=(a[:,0,0]-a[:,1,1]+1j*(a[:,0,1]+a[:,1,0]))/np.sqrt(2)
    ub=(b[:,0,0]-b[:,1,1]+1j*(b[:,0,1]+b[:,1,0]))/np.sqrt(2)
    ca=ka[:,0]+1j*ka[:,1];cb=kb[:,0]+1j*kb[:,1]
    z1=np.vdot(va,vb)+np.vdot(wa,wb)+2*np.vdot(ca,cb)/curvature_scale**2
    z2=np.vdot(ua,ub)
    aa,bb,cc,dd=z1.real,-z1.imag,z2.real,-z2.imag
    polynomial=np.array([-bb+2*dd,-2*aa+8*cc,-12*dd,-2*aa-8*cc,bb+2*dd])
    norm=np.max(abs(polynomial))
    candidates=[0.0,np.pi]
    if norm>1e-20:
        roots=np.roots(np.trim_zeros(polynomial/norm,'f'))
        candidates.extend(2*np.arctan(r.real) for r in roots if abs(r.imag)<1e-7*(1+abs(r.real)))
    scores=[aa*np.cos(p)+bb*np.sin(p)+cc*np.cos(2*p)+dd*np.sin(2*p) for p in candidates]
    return float(candidates[int(np.argmax(scores))])


def distance(qa,ka,qb,kb,curvature_scale,reduce=True):
    angle=best_angle(qa,ka,qb,kb,curvature_scale) if reduce else 0.0
    q,k=action(qb,kb,angle)
    angles=[]
    for a,b in zip(qa,q):
        a=a/np.linalg.norm(a);b=b/np.linalg.norm(b)
        angles.append(np.linalg.norm(model.logq(model.mul(np.r_[a[0],-a[1:]],b))))
    f=frames(qa)-frames(q)
    score=np.mean(np.sum(f*f,axis=(1,2))/2+np.sum((ka-k)**2,axis=1)/curvature_scale**2)
    return dict(angle=angle,rms=float(np.sqrt(score)),orientation_max=float(max(angles)),
                curvature_max=float(np.max(abs(ka-k))))


def check(output=None):
    rng=np.random.default_rng(260409488);n=10
    q=np.array([model.expq(v) for v in rng.normal(scale=.2,size=(n,3))])
    k=rng.normal(size=(n,2));angle=.731
    qt,kt=action(q,k,angle)
    beta,zeta,xi,chi,eta=model.parameters(lam=127,eta=1,chi=117)
    rold=model.geometry(model.rotate_increment(q,rng.normal(scale=.003,size=(n,3))))[1]
    rprev=model.geometry(model.rotate_increment(q,rng.normal(scale=.003,size=(n,3))))[1]
    kprev=k+rng.normal(scale=.05,size=k.shape)
    _,kp=action(q,kprev,angle)
    c,s=np.cos(angle),np.sin(angle);z=np.array([[c,-s,0.0],[s,c,0.0],[0.0,0.0,1.0]])
    r1=model.residual_q(q,k,rold,kprev,rprev,beta,zeta,xi,chi,eta,.006,2)[0]
    r2=model.residual_q(qt,kt,rold@z.T,kp,rprev@z.T,beta,zeta,xi,chi,eta,.006,2)[0]
    equivariance=float(np.max(abs(r2-r1@z.T))/np.max(abs(r1)))
    recovered=distance(q,k,qt,kt,chi)
    assert equivariance<1e-10 and recovered['rms']<1e-12
    # Sign choices must not affect orientation comparison.
    signed=qt.copy();signed[::2]*=-1
    assert distance(q,k,signed,kt,chi)['rms']<1e-12
    qb=np.array([model.expq(v) for v in rng.normal(scale=.8,size=(n,3))]);kb=rng.normal(scale=10,size=k.shape)
    optimum=distance(q,k,qb,kb,chi)
    def score(phi):
        qq,kk=action(qb,kb,phi);f=frames(q)-frames(qq)
        return float(np.mean(np.sum(f*f,axis=(1,2))/2+np.sum((k-kk)**2,axis=1)/chi**2))
    sampled=min(score(p) for p in np.linspace(-np.pi,np.pi,4097))
    assert optimum['rms']**2 <= sampled+1e-12
    result=dict(eta=1,residual_equivariance_relative_error=equivariance,
                recovered_rotation=recovered,quartic_minimum=optimum['rms']**2,minimum_4097_angle_grid=sampled,
                scope='Circular-section model only; full frames and curvatures compared, with quaternion signs identified.')
    if output is not None:
        Path(output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    print(json.dumps(check(Path(__file__).resolve().parent/'symmetry-check.json')))

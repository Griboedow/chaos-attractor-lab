"""Shen and Zhang, arXiv:2402.18363v1, literal smooth Eq. 6 / A1.

State: nm,nc,xo,y,xs. Time is the paper's dimensionless tau=t*omega_d0.
This equation extends the linear motor force-velocity branch past stall and
omits the absolute values of the preceding physical equations. Those are
explicit model limitations; no clipping or branch substitution is performed.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import lambertw

SOURCE='https://arxiv.org/abs/2402.18363v1'
NAMES=['Nm','Nc','omega','kon','koff','km','kc','ks','fs','fd','Fb','vu','Gamma/gamma','F']


def parameters(nc=4.,vu=1.5,ratio=1.,stiffness_factor=1.,omega=40/350,load=0.):
    # Rounded characteristic length and numerical force scale from Table I.
    # The force value 2.3033 is interpreted in pN, consistent with l0 and
    # the plotted equilibrium y=9.2436, not with its printed kBT/nm label.
    length=1.7887;force=2.3033
    km=.3*length/force;kc=.03144*length/force
    return np.array([100.,nc,omega,1/350,.1/350,km,kc,kc*stiffness_factor,
                     4.96/force,2.4/force,6.25/force,vu,ratio,load])


def field_jac(x,p):
    nm,nc,xo,y,xs=x
    Nm,Nc,w,on,off,km,kc,ks,fs,fd,Fb,vu,ratio,F=p
    xc=xo-xs;em=np.exp(km*y/fd);ec=np.exp(kc*xc/Fb)
    u=nm*km*y-nc*kc*xc-F
    f=np.array([w*(Nm-nm)-nm*em,on*(Nc-nc)-off*nc*ec,
                u,vu*(1-km*y/fs)-u,ratio*(nc*kc*xc-ks*xs)])
    j=np.zeros((5,5),dtype=np.result_type(x,p))
    j[0,0]=-w-em;j[0,3]=-nm*em*km/fd
    j[1,1]=-on-off*ec;j[1,2]=-off*nc*ec*kc/Fb;j[1,4]=-j[1,2]
    j[2]=[km*y,-kc*xc,-nc*kc,nm*km,nc*kc]
    j[3]=-j[2];j[3,3]-=vu*km/fs
    j[4]=[0.,ratio*kc*xc,ratio*nc*kc,0.,-ratio*(nc*kc+ks)]
    return f,j


def independent_field(x,p):
    """Force-balance form in clutch extension, then recover actin extension."""
    m,c,actin,motor,substrate=x
    Nm,Nc,w,on,off,km,kc,ks,fs,fd,Fb,vu,ratio,F=p
    motor_force=km*motor;clutch_force=kc*(actin-substrate)
    substrate_velocity=ratio*(c*clutch_force-ks*substrate)
    clutch_velocity=m*motor_force-(1+ratio)*c*clutch_force+ratio*ks*substrate-F
    actin_velocity=clutch_velocity+substrate_velocity
    return np.array([w*(Nm-m)-m*np.exp(motor_force/fd),
                     on*(Nc-c)-off*c*np.exp(clutch_force/Fb),actin_velocity,
                     vu*(1-motor_force/fs)-actin_velocity,substrate_velocity])


def equilibria(p):
    Nm,Nc,w,on,off,km,kc,ks,fs,fd,Fb,vu,ratio,F=p
    m=w*Nm/(w+np.exp(fs/fd));force=m*fs-F;a=force/Fb
    def g(c):return c+(off/on)*c*np.exp(a/c)-Nc
    if a>0:
        minimum=a/(1+lambertw(on/(off*np.e)).real)
        if minimum>=Nc or g(minimum)>0:return []
        roots=[brentq(g,max(a/600,1e-15),minimum),brentq(g,minimum,Nc)]
    else: roots=[brentq(g,1e-15,Nc)]
    result=[]
    for c in roots:
        substrate=force/ks;clutch=force/(c*kc)
        x=np.array([m,c,substrate+clutch,fs/km,substrate])
        f,j=field_jac(x,p)
        result.append(dict(state=x,eigenvalues=np.linalg.eigvals(j),residual=float(max(abs(f)))))
    return result

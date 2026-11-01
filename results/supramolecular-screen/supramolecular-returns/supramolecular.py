"""Sharko et al. 2022 simplified cooperative polymer model, Eqs5.6/5.12.

Also printed verbatim in Singh et al. arXiv2601.15662v1 Eqs1-4.
This is the approximate moment closure, not the full polymer size network.
Set k+=1. State(x,y,z)=(a1/kd2,d/kd2,m0/kd3), tau=kd3*t.
Parameters(alpha,beta,T,C)=(ka/kd3,kd2/kd3,ctot/kd2,C).
Polymer mass/kd2 is m=T-x-y; physical domain m>2*z/beta.
"""
import numpy as np
from numba import njit
from scipy.optimize import brentq

SOURCE='https://doi.org/10.1021/acs.chemrev.1c00958'
REFERENCE=np.array([.2,1000.,10.,4.6])


@njit(cache=True)
def field_jac(q,p):
    x,y,z=q;alpha,beta,total,C=p[:4];mass=total-x-y;gap=beta*mass-2*z
    f=np.array([alpha*y-2*x*z,-alpha*y+2*z,beta*mass-6*z-2*beta*C*z**3/gap**2-z*z])
    slope=1+4*beta*C*z**3/gap**3
    dz=-6-2*beta*C*(3*z*z/gap**2+4*z**3/gap**3)-2*z
    j=np.array([[-2*z,alpha,-2*x],[0.,-alpha,2.],[-beta*slope,-beta*slope,dz]])
    return f,j


def source_four_field(v,p):
    """Literal unscaled four-equation source, kd3=k+=1, kd2=beta."""
    d,a,mass,n=v;ka,kd2,total,C=p[:4]
    release=2*kd2*n;elongation=2*a*n
    return np.array([-ka*d+release,ka*d-elongation,-release+elongation,
        mass-4*n-2*kd2*C*n/(mass/n-2)**2-2*n-n*n])


def independent_field(q,p):
    alpha,beta,total,C=p[:4];x,y,z=q
    original=np.array([beta*y,beta*x,beta*(total-x-y),z])
    out=source_four_field(original,p)
    return np.array([out[1]/beta,out[0]/beta,out[3]])


def equilibrium(p):
    alpha,beta,total,C=p[:4]
    if total<=1:raise ValueError('Positive interior equilibrium requires T>1.')
    upper=beta*(total-1)/(2*beta/alpha+2)
    def fun(z):return field_jac(np.array([1.,2*z/alpha,z]),p)[0][2]
    z=brentq(fun,upper*1e-12,upper*(1-1e-10),xtol=1e-14,rtol=1e-14)
    state=np.array([1.,2*z/alpha,z]);f,j=field_jac(state,p)
    return dict(state=state,eigenvalues=np.linalg.eigvals(j),residual=float(np.max(abs(f))),length=float(beta*(total-state[0]-state[1])/z))

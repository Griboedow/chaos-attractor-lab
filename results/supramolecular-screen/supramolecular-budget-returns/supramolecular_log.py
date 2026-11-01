"""Positive polymer coordinates: (x, log(z), log(L-2)), L=beta*m/z.

Exact interior change of coordinates of supramolecular.py, not a new model.
Avoids subtracting two almost equal total monomer pools to obtain polymer mass.
Log concentrations can remain resolved when their exponentials underflow.
Finite-time tangent rates in these coordinates need not equal raw-state rates.
"""
import numpy as np
from numba import njit
from supramolecular import source_four_field


@njit(cache=True)
def field_jac(q,p):
    x,u,v=q;alpha,beta,total,C=p[:4]
    z=np.exp(u);excess=np.exp(v);length=2+excess;mass=z*length/beta
    loss=2*beta*C/excess**2
    du=excess-4-loss-z
    dv=(2*beta*(x-1)-length*du)/excess
    f=np.array([alpha*(total-x-mass)-2*x*z,du,dv])
    j=np.array([[-alpha-2*z,-alpha*mass-2*x*z,-alpha*z*excess/beta],
        [0.,-z,excess+2*loss],
        [2*beta/excess,length*z/excess,-du-length*(excess+2*loss)/excess-dv]])
    return f,j


def from_raw(q,p):
    x,y,z=q;mass=p[2]-x-y;length=p[1]*mass/z
    if min(x,y,z)<=0 or length<=2:raise ValueError('Start must be in the positive interior.')
    return np.array([x,np.log(z),np.log(length-2)])


def to_raw(q,p):
    q=np.asarray(q);x=q[...,0];z=np.exp(q[...,1]);mass=z*(2+np.exp(q[...,2]))/p[1]
    return np.stack([x,p[2]-x-mass,z],axis=-1)


def independent_field(q,p):
    x,u,v=q;alpha,beta,total,C=p[:4]
    n=np.exp(u);length=2+np.exp(v);mass=n*length
    a=beta*x;d=beta*total-a-mass
    source=source_four_field(np.array([d,a,mass,n]),p)
    du=source[3]/n
    return np.array([source[1]/beta,du,length/(length-2)*(source[2]/mass-du)])

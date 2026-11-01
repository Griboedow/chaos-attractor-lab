"""Separately implemented compiled field and fixed-step local-return controls."""
import numpy as np
from numba import njit
from icosa_poincare import c,u,s,g,x0,n,B,inv

@njit
def rhs(t,x):
    f=np.zeros(3);J=-g*np.eye(3);ss=s*s
    for j in range(60):
        a=x[0]-c[j,0];b=x[1]-c[j,1];e=x[2]-c[j,2];w=np.exp(-.5*(a*a+b*b+e*e)/ss)
        p=u[j,1]*e-u[j,2]*b;q=u[j,2]*a-u[j,0]*e;r=u[j,0]*b-u[j,1]*a
        f[0]+=w*p;f[1]+=w*q;f[2]+=w*r
        J[0,0]-=w*p*a/ss;J[0,1]+=w*(-u[j,2]-p*b/ss);J[0,2]+=w*(u[j,1]-p*e/ss)
        J[1,0]+=w*(u[j,2]-q*a/ss);J[1,1]-=w*q*b/ss;J[1,2]+=w*(-u[j,0]-q*e/ss)
        J[2,0]+=w*(-u[j,1]-r*a/ss);J[2,1]+=w*(u[j,0]-r*b/ss);J[2,2]-=w*r*e/ss
    f-=g*x[:3]
    if len(x)==3:return f
    out=np.empty(9);out[:3]=f;out[3:]=(J@x[3:].reshape(3,2)).ravel();return out
@njit
def step(x,h):
    a=rhs(0,x);b=rhs(0,x+.5*h*a);e=rhs(0,x+.5*h*b);z=rhs(0,x+h*e);return x+h/6*(a+2*b+2*e+z)
@njit
def rk_return(uv,h,horizon=25.):
    x=x0+B@uv;old=(x-x0)@n;t=0.
    for _ in range(int(np.ceil(horizon/h))):
        y=step(x,h);new=(y-x0)@n
        if old<0 and new>=0 and t>1.:
            lo=0.;hi=h
            for _ in range(32):
                mid=.5*(lo+hi);z=step(x,mid)
                if (z-x0)@n<0:lo=mid
                else:hi=mid
            z=step(x,.5*(lo+hi));q=inv@(z-x0);speed=n@rhs(0,z)
            if abs(q[0])<.15 and abs(q[1])<.01 and speed>.05:return q,t+.5*(lo+hi),speed
        x=y;old=new;t+=h
    return np.full(2,np.nan),np.nan,np.nan

"""Steady localized Gaussian vortices, linear confinement and constant drift.

p = gamma, omega[3], sigma[3], center_distance, drift[3].
Centers are center_distance*e_n; vortex axes are e_(n+1 mod 3).
"""
import numpy as np
from numba import njit, prange
from twist_map import orthonormalize


@njit
def field_jacobian(z,p):
    f=-p[0]*z+p[8:11]
    jac=-p[0]*np.eye(3)
    for n in range(3):
        r=z.copy();r[n]-=p[7]
        axis=(n+1)%3;i=(axis+1)%3;j=(axis+2)%3
        sigma2=p[4+n]**2
        g=p[1+n]*np.exp(-.5*np.dot(r,r)/sigma2)
        cross=np.zeros(3);cross[i]=-r[j];cross[j]=r[i]
        f+=g*cross
        for a in range(3):
            for b in range(3):
                jac[a,b]-=g*cross[a]*r[b]/sigma2
        jac[i,j]-=g;jac[j,i]+=g
    return f,jac


@njit
def advance(z,q,p,h):
    a,j=field_jacobian(z,p);av=j@q
    b,j=field_jacobian(z+h*.5*a,p);bv=j@(q+h*.5*av)
    c,j=field_jacobian(z+h*.5*b,p);cv=j@(q+h*.5*bv)
    d,j=field_jacobian(z+h*c,p);dv=j@(q+h*cv)
    return z+h/6*(a+2*b+2*c+d),q+h/6*(av+2*bv+2*cv+dv)


@njit(parallel=True)
def integrate(parameters,starts,stepsizes,warm_time,observed_time,window_time,keep,stride):
    count=len(parameters);windows=np.zeros((count,int(observed_time/window_time),3))
    tails=np.zeros((count,keep,3));ends=np.zeros((count,3));means=np.zeros_like(windows)
    max_radius=np.zeros(count)
    for lane in prange(count):
        h=stepsizes[lane];warm=int(round(warm_time/h));observed=int(round(observed_time/h));window=int(round(window_time/h))
        z=starts[lane].copy();q=np.eye(3);acc=np.zeros(3);average=np.zeros(3)
        for t in range(warm+observed):
            z,q=advance(z,q,parameters[lane],h);q,logs=orthonormalize(q)
            if t>=warm:
                acc+=logs;average+=z
                max_radius[lane]=max(max_radius[lane],np.sqrt(np.dot(z,z)))
                if (t-warm+1)%window==0:
                    ii=(t-warm+1)//window-1;windows[lane,ii]=acc/(window*h);means[lane,ii]=average/window
                    acc[:]=0;average[:]=0
                offset=t-(warm+observed-keep*stride)
                if offset>=0 and offset%stride==0:tails[lane,offset//stride]=z
        ends[lane]=z
    return windows,tails,ends,means,max_radius


@njit
def trajectory(p,start,h,warm,observed,stride):
    z=start.copy();q=np.zeros((3,1))
    for _ in range(warm):z,q=advance(z,q,p,h)
    points=np.empty((observed//stride,3))
    for i in range(observed):
        z,q=advance(z,q,p,h)
        if (i+1)%stride==0:points[(i+1)//stride-1]=z
    return points

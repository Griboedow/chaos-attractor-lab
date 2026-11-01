"""Published multiplicatively overrelaxed Sinkhorn in anchored log coordinates."""
import numpy as np
from numba import njit

@njit(cache=True)
def lse(a):
    peak=np.max(a)
    return peak+np.log(np.sum(np.exp(a-peak)))

@njit(cache=True)
def targets(logk):
    m,n=logk.shape
    a=np.empty(m);b=np.empty(n)
    for i in range(m):a[i]=lse(logk[i])
    for j in range(n):b[j]=lse(logk[:,j])
    return a,b

@njit(cache=True)
def step(z,v,logk,loga,logb,omega):
    m,n=logk.shape;split=m-1
    x=np.zeros(m);y=np.zeros(n);vx=np.zeros(m);vy=np.zeros(n)
    x[1:]=z[:split];y[1:]=z[split:]
    vx[1:]=v[:split];vy[1:]=v[split:]
    f=np.empty(m);df=np.empty(m)
    for i in range(m):
        logits=logk[i]+y;normal=lse(logits)
        f[i]=loga[i]-normal
        df[i]=-np.sum(np.exp(logits-normal)*vy)
    a=1-omega
    for i in range(1,m):
        x[i]=a*x[i]+omega*(f[i]-f[0])
        vx[i]=a*vx[i]+omega*(df[i]-df[0])
    g=np.empty(n);dg=np.empty(n)
    for j in range(n):
        logits=logk[:,j]+x;normal=lse(logits)
        g[j]=logb[j]-normal
        dg[j]=-np.sum(np.exp(logits-normal)*vx)
    for j in range(1,n):
        y[j]=a*y[j]+omega*(g[j]-g[0])
        vy[j]=a*vy[j]+omega*(dg[j]-dg[0])
    return np.concatenate((x[1:],y[1:])),np.concatenate((vx[1:],vy[1:]))

def jacobian(z,logk,loga,logb,omega):
    return np.column_stack([step(z,v,logk,loga,logb,omega)[1] for v in np.eye(len(z))])

def bounds(logk,loga,logb,omega):
    rows=logk[1:]-logk[0]
    cols=(logk[:,1:]-logk[:,0,None]).T
    lo=np.r_[loga[1:]-loga[0]-rows.max(axis=1),logb[1:]-logb[0]-cols.max(axis=1)]
    hi=np.r_[loga[1:]-loga[0]-rows.min(axis=1),logb[1:]-logb[0]-cols.min(axis=1)]
    forcing=np.maximum(abs(lo),abs(hi))
    return omega*forcing/(1-abs(1-omega))

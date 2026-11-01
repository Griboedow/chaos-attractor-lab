"""Unequal fixed block relaxation; an explicit extension of the single-rate map."""
import numpy as np
from numba import njit
from model import lse,targets

@njit(cache=True)
def step(z,v,L,la,lb,wx,wy):
    m,n=L.shape;p=m-1
    x=np.zeros(m);y=np.zeros(n);vx=np.zeros(m);vy=np.zeros(n)
    x[1:]=z[:p];y[1:]=z[p:];vx[1:]=v[:p];vy[1:]=v[p:]
    f=np.empty(m);df=np.empty(m)
    for i in range(m):
        logits=L[i]+y;s=lse(logits);f[i]=la[i]-s;df[i]=-np.sum(np.exp(logits-s)*vy)
    for i in range(1,m):
        x[i]=(1-wx)*x[i]+wx*(f[i]-f[0]);vx[i]=(1-wx)*vx[i]+wx*(df[i]-df[0])
    g=np.empty(n);dg=np.empty(n)
    for j in range(n):
        logits=L[:,j]+x;s=lse(logits);g[j]=lb[j]-s;dg[j]=-np.sum(np.exp(logits-s)*vx)
    for j in range(1,n):
        y[j]=(1-wy)*y[j]+wy*(g[j]-g[0]);vy[j]=(1-wy)*vy[j]+wy*(dg[j]-dg[0])
    return np.concatenate((x[1:],y[1:])),np.concatenate((vx[1:],vy[1:]))

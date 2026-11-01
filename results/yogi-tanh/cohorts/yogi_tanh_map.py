"""Normalized dense TensorFlow Addons Yogi with its existing tanh activation.

For f(theta)=theta**2/2, alpha=sqrt(k/10), epsilon=alpha*e,
(x,m,v)=(theta/alpha, physical_m/alpha, physical_v/alpha**2).
The library's tanh(10*(g*g-v)) becomes tanh(k*(x*x-v)).
This is a source transcription, not an upstream-library execution.
"""
import numpy as np
from numba import njit, prange


@njit
def step(z,p,correction=1.):
    x,m,v=z; b,c,e,k=p
    T=np.tanh(k*(v-x*x)); M=b*m+(1-b)*x
    V=v-(1-c)*T*x*x
    return np.array([x-correction*M/(np.sqrt(V)+e),M,V])


@njit
def step_tangent(z,q,p,correction=1.):
    x,m,v=z; b,c,e,k=p
    T=np.tanh(k*(v-x*x)); slope=1-T*T
    M=b*m+(1-b)*x; V=v-(1-c)*T*x*x
    Vx=-2*(1-c)*x*(T-k*x*x*slope)
    Vv=1-(1-c)*k*x*x*slope
    root=np.sqrt(V); d=root+e; result=np.empty_like(q)
    for j in range(q.shape[1]):
        dm=(1-b)*q[0,j]+b*q[1,j]; dv=Vx*q[0,j]+Vv*q[2,j]
        result[0,j]=q[0,j]-correction*dm/d+correction*M*dv/(2*root*d*d)
        result[1,j]=dm; result[2,j]=dv
    return np.array([x-correction*M/d,M,V]),result,b*Vv


@njit(parallel=True)
def scout(parameters,starts,bias,warm,observed,window,keep):
    n=len(parameters); windows=np.full((n,observed//window),np.nan)
    tails=np.full((n,keep,3),np.nan); ends=starts.copy()
    status=np.zeros(n,dtype=np.int64); extrema=np.full((n,2,3),np.nan)
    for i in prange(n):
        b,c,e,k=parameters[i];x,m,v=starts[i];qx,qm,qv=.3,.7,.4
        norm=np.sqrt(qx*qx+qm*qm+qv*qv);qx/=norm;qm/=norm;qv/=norm
        pb=1.;pc=1.;acc=0.;lo=np.full(3,np.inf);hi=np.full(3,-np.inf)
        for t in range(warm+observed):
            pb=pb*b if pb>1e-300 else 0.;pc=pc*c if pc>1e-300 else 0.
            factor=np.sqrt(1-pc)/(1-pb) if bias[i] else 1.
            T=np.tanh(k*(v-x*x));slope=1-T*T
            M=b*m+(1-b)*x;V=v-(1-c)*T*x*x
            if V<=0:status[i]=-(t+1);break
            root=np.sqrt(V);d=root+e
            Vx=-2*(1-c)*x*(T-k*x*x*slope);Vv=1-(1-c)*k*x*x*slope
            dm=(1-b)*qx+b*qm;dv=Vx*qx+Vv*qv
            dx=qx-factor*dm/d+factor*M*dv/(2*root*d*d)
            X=x-factor*M/d;norm=np.sqrt(dx*dx+dm*dm+dv*dv)
            if not np.isfinite(X+M+V+norm) or max(abs(X),abs(M),abs(V))>1e8 or norm<1e-280:
                status[i]=t+1;break
            x,m,v=X,M,V;qx,qm,qv=dx/norm,dm/norm,dv/norm
            if t>=warm:
                acc+=np.log(norm);j=t-warm+1
                for a,value in enumerate((x,m,v)):lo[a]=min(lo[a],value);hi[a]=max(hi[a],value)
                if j%window==0:windows[i,j//window-1]=acc/window;acc=0.
                if j>observed-keep:tails[i,j-observed+keep-1]=np.array([x,m,v])
        ends[i]=np.array([x,m,v]);extrema[i,0]=lo;extrema[i,1]=hi
    return windows,tails,ends,status,extrema


@njit(parallel=True)
def spectra(parameters,starts,bias,warm,observed,window,keep,Q):
    n=len(parameters);windows=np.full((n,observed//window,3),np.nan)
    det_windows=np.full((n,observed//window),np.nan);tails=np.full((n,keep,3),np.nan)
    ends=starts.copy();status=np.zeros(n,dtype=np.int64);negative=np.zeros(n,dtype=np.int64)
    for i in prange(n):
        p=parameters[i];z=starts[i].copy();q=Q.copy();pb=1.;pc=1.;acc=np.zeros(3);detacc=0.
        for t in range(warm+observed):
            pb=pb*p[0] if pb>1e-300 else 0.;pc=pc*p[1] if pc>1e-300 else 0.
            factor=np.sqrt(1-pc)/(1-pb) if bias[i] else 1.
            z,q,det=step_tangent(z,q,p,factor)
            if not np.isfinite(z).all() or not np.isfinite(q).all() or det==0:status[i]=t+1;break
            logs=np.empty(3)
            for j in range(3):
                for _ in range(2):
                    for h in range(j):q[:,j]-=np.dot(q[:,j],q[:,h])*q[:,h]
                norm=np.linalg.norm(q[:,j])
                if norm==0 or not np.isfinite(norm):status[i]=-(t+1);break
                q[:,j]/=norm;logs[j]=np.log(norm)
            if status[i]:break
            if t>=warm:
                acc+=logs;detacc+=np.log(abs(det));negative[i]+=int(det<0);j=t-warm+1
                if j%window==0:windows[i,j//window-1]=acc/window;det_windows[i,j//window-1]=detacc/window;acc[:]=0.;detacc=0.
                if j>observed-keep:tails[i,j-observed+keep-1]=z
        ends[i]=z
    return windows,det_windows,tails,ends,status,negative


@njit
def trajectory(p,start,warm,observed):
    z=start.copy()
    for _ in range(warm):z=step(z,p)
    states=np.empty((observed,3))
    for i in range(observed):z=step(z,p);states[i]=z
    return states

"""Scalar-quadratic Adamax with the PyTorch 2.15 epsilon placement.

For f(theta)=theta^2/2, x=theta/alpha,m=M/alpha,u=U/alpha,e=epsilon/alpha.
State order (x,m,u), parameters (beta1,beta2,e). The autonomous limit has
correction=1. Original bias-corrected updates use correction=1-beta1**t.
"""
import numpy as np
from numba import njit,prange


@njit
def step_tangent(z,q,p,correction=1.):
    x,m,u=z;b,c,e=p
    mn=b*m+(1-b)*x;reset=abs(x)+e>=c*u
    un=abs(x)+e if reset else c*u
    result=np.empty_like(q)
    for j in range(q.shape[1]):
        dm=b*q[1,j]+(1-b)*q[0,j]
        du=(1. if x>0 else -1. if x<0 else 0.)*q[0,j] if reset else c*q[2,j]
        result[0,j]=q[0,j]-dm/(correction*un)+mn*du/(correction*un*un)
        result[1,j]=dm;result[2,j]=du
    return np.array([x-mn/(correction*un),mn,un]),result,reset


@njit
def step(z,p,correction=1.):
    x,m,u=z;b,c,e=p;mn=b*m+(1-b)*x;un=max(c*u,abs(x)+e)
    return np.array([x-mn/(correction*un),mn,un])


@njit(parallel=True)
def scout(parameters,starts,warm,observed,window,keep):
    count=len(parameters);windows=np.full((count,observed//window),np.nan)
    tails=np.full((count,keep,3),np.nan);ends=np.zeros((count,3));status=np.zeros(count,dtype=np.int32)
    extrema=np.zeros((count,2,3));resets=np.zeros(count,dtype=np.int64);margins=np.full(count,np.inf)
    for lane in prange(count):
        b,c,e=parameters[lane];x,m,u=starts[lane];vx,vm,vu=1.,0.,0.;power=1.;acc=0.
        lo=np.full(3,np.inf);hi=np.full(3,-np.inf)
        for t in range(warm+observed):
            power*=b;correction=1-power if t<warm else 1.
            mn=b*m+(1-b)*x;dm=b*vm+(1-b)*vx
            margin=abs(x)+e-c*u;reset=margin>=0
            if reset:
                un=abs(x)+e;du=(1. if x>0 else -1. if x<0 else 0.)*vx
            else:un=c*u;du=c*vu
            xn=x-mn/(correction*un);dx=vx-dm/(correction*un)+mn*du/(correction*un*un)
            norm=np.sqrt(dx*dx+dm*dm+du*du)
            if not np.isfinite(xn+mn+un+norm) or max(abs(xn),abs(mn),abs(un))>1e6 or norm<1e-280:
                status[lane]=t+1;break
            x,m,u=xn,mn,un;vx,vm,vu=dx/norm,dm/norm,du/norm
            if t>=warm:
                acc+=np.log(norm);resets[lane]+=int(reset);margins[lane]=min(margins[lane],abs(margin)/(1+un))
                for j,value in enumerate((x,m,u)):lo[j]=min(lo[j],value);hi[j]=max(hi[j],value)
                n=t-warm+1
                if n%window==0:windows[lane,n//window-1]=acc/window;acc=0.
                if n>observed-keep:tails[lane,n-(observed-keep)-1]=np.array([x,m,u])
        ends[lane]=np.array([x,m,u]);extrema[lane,0]=lo;extrema[lane,1]=hi
    return windows,tails,ends,status,extrema,resets,margins


@njit
def trajectory(p,start,warm,observed):
    z=start.copy()
    for _ in range(warm):z=step(z,p)
    data=np.empty((observed,3))
    for i in range(observed):z=step(z,p);data[i]=z
    return data


@njit(parallel=True)
def spectra(parameters,starts,warm,observed,window,keep,bias_warm=False):
    """Two surviving tangent exponents. A reset makes the third exactly -infinity."""
    count=len(parameters);windows=np.zeros((count,observed//window,2));tails=np.zeros((count,keep,3))
    ends=np.zeros((count,3));reset_counts=np.zeros(count,dtype=np.int64);minimum_margins=np.full(count,np.inf)
    for lane in prange(count):
        z=starts[lane].copy();q=np.eye(3)[:,:2].copy();acc=np.zeros(2);p=parameters[lane];power=1.
        for t in range(warm+observed):
            power*=p[0];correction=1-power if bias_warm and t<warm else 1.
            margin=abs(abs(z[0])+p[2]-p[1]*z[2])/(1+z[2]);z,q,reset=step_tangent(z,q,p,correction)
            logs=np.empty(2)
            for j in range(2):
                for _ in range(2):
                    for k in range(j):q[:,j]-=np.dot(q[:,j],q[:,k])*q[:,k]
                norm=np.linalg.norm(q[:,j]);q[:,j]/=norm;logs[j]=np.log(norm)
            if t>=warm:
                acc+=logs;reset_counts[lane]+=int(reset);minimum_margins[lane]=min(minimum_margins[lane],margin)
                n=t-warm+1
                if n%window==0:windows[lane,n//window-1]=acc/window;acc[:]=0.
                if n>observed-keep:tails[lane,n-(observed-keep)-1]=z
        ends[lane]=z
    return windows,tails,ends,reset_counts,minimum_margins

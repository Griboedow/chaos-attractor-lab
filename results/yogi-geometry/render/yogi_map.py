"""Scalar-quadratic Yogi: the published autonomous map, without bias correction.

For f(theta)=theta**2/2, x=theta/alpha, m=M/alpha, v=V/alpha**2,
and e=epsilon/alpha. This linear normalization is not an artistic deformation.
The map is discontinuous at v=x*x; branch derivatives omit that boundary.
"""
import numpy as np
from numba import njit, prange


@njit
def step_tangent(z, q, p, correction_m=1., correction_v=1.):
    x,m,v=z; b,c,e=p
    sign=1. if v>x*x else -1. if v<x*x else 0.
    mn=b*m+(1-b)*x; vn=v-(1-c)*sign*x*x
    root=np.sqrt(vn/correction_v); denominator=root+e
    result=np.empty_like(q)
    for j in range(q.shape[1]):
        dm=b*q[1,j]+(1-b)*q[0,j]
        dv=q[2,j]-2*(1-c)*sign*x*q[0,j]
        result[0,j]=q[0,j]-dm/(correction_m*denominator)+mn*dv/(2*correction_m*correction_v*root*denominator**2)
        result[1,j]=dm; result[2,j]=dv
    return np.array([x-mn/(correction_m*denominator),mn,vn]),result


@njit
def step(z,p,correction_m=1.,correction_v=1.):
    x,m,v=z; b,c,e=p
    sign=1. if v>x*x else -1. if v<x*x else 0.
    mn=b*m+(1-b)*x; vn=v-(1-c)*sign*x*x
    return np.array([x-mn/(correction_m*(np.sqrt(vn/correction_v)+e)),mn,vn])


@njit(parallel=True)
def scout(parameters,starts,warm,observed,window,keep):
    count=len(parameters); windows=np.full((count,observed//window),np.nan)
    tails=np.full((count,keep,3),np.nan); ends=np.zeros((count,3)); status=np.zeros(count,dtype=np.int32)
    extrema=np.zeros((count,2,3)); margins=np.full(count,np.inf); changes=np.zeros(count,dtype=np.int64)
    for lane in prange(count):
        b,c,e=parameters[lane]; x,m,v=starts[lane]; vx,vm,vv=1.,0.,0.; acc=0.; oldsign=0.
        lo=np.full(3,np.inf); hi=np.full(3,-np.inf)
        for t in range(warm+observed):
            margin=v-x*x; sign=1. if margin>0 else -1. if margin<0 else 0.
            mn=b*m+(1-b)*x; vn=v-(1-c)*sign*x*x
            if vn<=0:
                status[lane]=-(t+1); break
            dm=b*vm+(1-b)*vx; dv=vv-2*(1-c)*sign*x*vx
            root=np.sqrt(vn); den=root+e
            xn=x-mn/den; dx=vx-dm/den+mn*dv/(2*root*den*den)
            norm=np.sqrt(dx*dx+dm*dm+dv*dv)
            if not np.isfinite(xn+mn+vn+norm) or max(abs(xn),abs(mn),abs(vn))>1e8 or norm<1e-280:
                status[lane]=t+1; break
            x,m,v=xn,mn,vn; vx,vm,vv=dx/norm,dm/norm,dv/norm
            if t>=warm:
                acc+=np.log(norm); margins[lane]=min(margins[lane],abs(margin)/(1+v))
                changes[lane]+=int(sign!=oldsign)
                for j,value in enumerate((x,m,v)):lo[j]=min(lo[j],value);hi[j]=max(hi[j],value)
                n=t-warm+1
                if n%window==0:windows[lane,n//window-1]=acc/window;acc=0.
                if n>observed-keep:tails[lane,n-(observed-keep)-1]=np.array([x,m,v])
            oldsign=sign
        ends[lane]=np.array([x,m,v]);extrema[lane,0]=lo;extrema[lane,1]=hi
    return windows,tails,ends,status,extrema,margins,changes


@njit
def trajectory(p,start,warm,observed):
    z=start.copy()
    for _ in range(warm):z=step(z,p)
    data=np.empty((observed,3))
    for i in range(observed):z=step(z,p);data[i]=z
    return data


@njit(parallel=True)
def spectra(parameters,starts,warm,observed,window,keep,bias_warm=False):
    count=len(parameters);windows=np.zeros((count,observed//window,3));tails=np.zeros((count,keep,3))
    ends=np.zeros((count,3));minimum_margins=np.full(count,np.inf)
    for lane in prange(count):
        z=starts[lane].copy();q=np.eye(3);acc=np.zeros(3);p=parameters[lane];pb=1.;pc=1.
        for t in range(warm+observed):
            pb*=p[0];pc*=p[1]
            cm=1-pb if bias_warm and t<warm else 1.;cv=1-pc if bias_warm and t<warm else 1.
            margin=abs(z[2]-z[0]**2)/(1+z[2]);z,q=step_tangent(z,q,p,cm,cv)
            logs=np.empty(3)
            for j in range(3):
                for _ in range(2):
                    for k in range(j):q[:,j]-=np.dot(q[:,j],q[:,k])*q[:,k]
                norm=np.linalg.norm(q[:,j]);q[:,j]/=norm;logs[j]=np.log(norm)
            if t>=warm:
                acc+=logs;minimum_margins[lane]=min(minimum_margins[lane],margin)
                n=t-warm+1
                if n%window==0:windows[lane,n//window-1]=acc/window;acc[:]=0.
                if n>observed-keep:tails[lane,n-(observed-keep)-1]=z
        ends[lane]=z
    return windows,tails,ends,minimum_margins

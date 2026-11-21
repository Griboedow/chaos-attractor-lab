"""Native-coordinate binary64 map, analytic Jacobian, and full three-vector QR."""
import math
import numpy as np
from numba import njit


@njit(cache=True)
def step(x,raw):
    beta,inv,tube,nf=raw[0]
    n=int(nf);w=np.empty(n);m=np.zeros(3);scale=raw[2:5,3]
    for k in range(n):
        w[k]=-inv*np.sum((x-raw[2+k,:3])**2)
    w=np.exp(w-w.max());w/=w.sum()
    for k in range(n):
        m+=w[k]*raw[2+k,:3]
    cov=np.zeros((3,3))
    for k in range(n):
        for a in range(3):
            for b in range(3):
                cov[a,b]+=2*inv*w[k]*(raw[2+k,a]-m[a])*(raw[2+k,b]-m[b])
    v=(2*m-x)/scale
    warp=raw[5,3]
    if warp!=0:v[2]-=warp*(v[0]**2-v[1]**2)
    rho=math.sqrt(v[0]**2+v[1]**2)
    length=math.sqrt((rho-1)**2+v[2]**2)
    cs=v[0]/rho;sn=v[1]/rho;nr=(rho-1)/length;nz=v[2]/length
    p=np.array([(1+tube*nr)*cs,(1+tube*nr)*sn,tube*nz])*scale
    az=np.array([-sn,cs,0.]);mer=np.array([nz*cs,nz*sn,-nr])
    dp=np.empty((3,3))
    for a in range(3):
        for b in range(3):
            dp[a,b]=(az[a]*az[b]*(1+tube*nr)/rho+mer[a]*mer[b]*tube/length)*scale[a]/scale[b]
    if warp!=0:
        # Compose the inverse shear at the reflected point, followed by the
        # forward shear at the unwarped toroidal projection.
        u=p/scale
        inverse=np.eye(3);inverse[2,0]=-2*scale[2]*warp*v[0]/scale[0];inverse[2,1]=2*scale[2]*warp*v[1]/scale[1]
        forward=np.eye(3);forward[2,0]=2*scale[2]*warp*u[0]/scale[0];forward[2,1]=-2*scale[2]*warp*u[1]/scale[1]
        dp=forward@dp@inverse
        p[2]+=scale[2]*warp*(u[0]**2-u[1]**2)
    ident=np.eye(3)
    jac=beta*(ident+dp@(2*cov-ident))+(1-2*beta)*cov
    return beta*(x+p)+(1-2*beta)*m,jac


@njit(cache=True)
def orbit(raw,state,burn,window,windows,saved):
    x=state.copy();q=np.eye(3);rates=np.zeros((windows,3));cloud=np.empty((saved,3))
    total=burn+window*windows
    for i in range(total):
        x,j=step(x,raw);v=j@q
        for k in range(3):
            for repeat in range(2):
                for l in range(k):
                    gain=v[0,k]*v[0,l]+v[1,k]*v[1,l]+v[2,k]*v[2,l]
                    for a in range(3):v[a,k]-=gain*v[a,l]
            norm=math.sqrt(v[0,k]**2+v[1,k]**2+v[2,k]**2)
            if norm<1e-290:raise ValueError('Degenerate tangent vector')
            for a in range(3):v[a,k]/=norm
            if i>=burn:rates[(i-burn)//window,k]+=math.log(norm)/window
        q=v
        if i>=total-saved:cloud[i-total+saved]=x
    return cloud,rates


def independent(x,raw):
    b=raw[2:2+int(raw[0,3]),:3];delta=raw[0,1]*np.sum((b-x)**2,axis=1)
    w=np.exp(delta.min()-delta);m=w@b/w.sum();a=raw[2:5,3];v=2*m-x;u=v/a
    warp=raw[5,3]
    if warp!=0:u[2]-=warp*(u[0]**2-u[1]**2)
    center=np.array([u[0],u[1],0.])/np.linalg.norm(u[:2])
    standard=center+raw[0,2]*(u-center)/np.linalg.norm(u-center)
    if warp!=0:standard[2]+=warp*(standard[0]**2-standard[1]**2)
    point=a*standard
    return raw[0,0]*(2*point-v+x)/2+(1-raw[0,0])*m

"""Independent binary64 interval replay with directed polynomial trigonometry.

No system sine or cosine routine contributes to a verified enclosure. Machin's
identity supplies rational bounds for pi/2. The reduced-angle Taylor expansion
uses interval coefficients, interval Horner arithmetic and an explicit remainder.
"""
import argparse
import base64
import json
import math
import zlib
from fractions import Fraction
from pathlib import Path
from search import CACHE
import numpy as np
from numba import njit


def atan_bounds(q,n):
    s=sum(Fraction((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(n))
    t=s+Fraction((-1)**n,(2*n+1)*q**(2*n+1))
    return min(s,t),max(s,t)


def enclosure(x):
    f=float(x)
    if Fraction.from_float(f)==x:
        return (f,f)
    return (np.nextafter(f,-np.inf),np.nextafter(f,np.inf))


def trig_constants():
    a,b=atan_bounds(5,22);c,d=atan_bounds(239,9)
    # pi/2 = 8 atan(1/5) - 2 atan(1/239).
    pi2=(enclosure(8*a-2*d)[0],enclosure(8*b-2*c)[1])
    coeff=np.array([[enclosure(Fraction((-1)**k,math.factorial(2*k+j))) for k in range(10)] for j in (0,1)])
    # |r|<=1.1. Taylor's theorem with all real sine/cosine derivatives <=1.
    bound=Fraction(11,10)**20/math.factorial(20)
    bound=max(bound,Fraction(11,10)**19/math.factorial(19))
    remainder=enclosure(bound)[1]
    return pi2,coeff,remainder


@njit(inline="always")
def plus(a,b):
    lo,hi=a[0]+b[0],a[1]+b[1]
    if np.isnan(lo) or np.isnan(hi):
        return -np.inf,np.inf
    return np.nextafter(lo,-np.inf),np.nextafter(hi,np.inf)


@njit(inline="always")
def minus(a,b):
    lo,hi=a[0]-b[1],a[1]-b[0]
    if np.isnan(lo) or np.isnan(hi):
        return -np.inf,np.inf
    return np.nextafter(lo,-np.inf),np.nextafter(hi,np.inf)


@njit(inline="always")
def times(a,b):
    aa=a[0]*b[0];ab=a[0]*b[1];ba=a[1]*b[0];bb=a[1]*b[1]
    if np.isnan(aa) or np.isnan(ab) or np.isnan(ba) or np.isnan(bb):
        return -np.inf,np.inf
    return np.nextafter(min(aa,ab,ba,bb),-np.inf),np.nextafter(max(aa,ab,ba,bb),np.inf)


@njit(inline="always")
def maximum_abs(a):
    return max(abs(a[0]),abs(a[1]))


@njit(inline="always")
def overlap(a,b):
    l,h=max(a[0],b[0]),min(a[1],b[1])
    if l>h:
        raise ValueError("Inconsistent interval enclosures")
    return l,h


@njit(cache=True)
def sincos(t,pi2,cs,eps):
    x=t[0]/2+t[1]/2
    if not np.isfinite(x) or abs(x)>1e8:
        return (-1.,1.),(-1.,1.)
    radius=np.nextafter(max(abs(t[0]-x),abs(t[1]-x)),np.inf)
    if radius>=2:
        return (-1.,1.),(-1.,1.)
    n=round(x/pi2[0])
    r=minus((x,x),times((float(n),float(n)),pi2))
    if maximum_abs(r)>1.:
        return (-1.,1.),(-1.,1.)
    sq=times(r,r)
    co=(cs[0,9,0],cs[0,9,1]);si=(cs[1,9,0],cs[1,9,1])
    for k in range(8,-1,-1):
        co=plus(times(co,sq),cs[0,k]);si=plus(times(si,sq),cs[1,k])
    si=plus(times(si,r),(-eps,eps));co=plus(co,(-eps,eps))
    quadrant=n%4
    if quadrant==1:
        si,co=co,(-si[1],-si[0])
    elif quadrant==2:
        si,co=(-si[1],-si[0]),(-co[1],-co[0])
    elif quadrant==3:
        si,co=(-co[1],-co[0]),si
    # First-order Taylor around x, together with the global Lipschitz bound.
    sr=min(radius,times((radius,radius),plus((maximum_abs(co),maximum_abs(co)),(0.,np.nextafter(radius/2,np.inf))))[1])
    cr=min(radius,times((radius,radius),plus((maximum_abs(si),maximum_abs(si)),(0.,np.nextafter(radius/2,np.inf))))[1])
    s=plus(si,(-sr,sr));c=plus(co,(-cr,cr))
    return overlap(s,(-1.,1.)),overlap(c,(-1.,1.))


@njit(cache=True)
def rotate(state,p,pi2,cs,eps,tight=False):
    x=state.copy()
    columns=x.shape[1]
    for i in range(3):
        for j in range(columns):
            x[i,j]=times(x[i,j],p[0])
    x[0,0]=plus(x[0,0],minus((1.,1.),p[0]))
    for i in range(3):
        j,k=(i+1)%3,(i+2)%3
        w=x[k,0];sq=times(w,w)
        angle=plus(p[1+i],times(w,plus(p[4+i],times(sq,p[7+i]))))
        si,co=sincos(angle,pi2,cs,eps)
        u=x[i,0].copy();v=x[j,0].copy()
        nu=minus(times(u,co),times(v,si));nv=plus(times(u,si),times(v,co))
        deriv=plus(p[4+i],times((3.,3.),times(sq,p[7+i])))
        for a in range(1,columns):
            du=x[i,a].copy();dv=x[j,a].copy()
            da=times(deriv,x[k,a])
            # Product-rule expansion of a rotation with a variable angle.
            dsi=times(co,da);dco=times((-si[1],-si[0]),da)
            au=minus(plus(times(du,co),times(u,dco)),plus(times(dv,si),times(v,dsi)))
            av=plus(plus(times(du,si),times(u,dsi)),plus(times(dv,co),times(v,dco)))
            if tight:
                bu=minus(minus(times(co,du),times(si,dv)),times(nv,da))
                bv=plus(plus(times(si,du),times(co,dv)),times(nu,da))
                x[i,a]=overlap(au,bu);x[j,a]=overlap(av,bv)
            else:
                x[i,a]=au;x[j,a]=av
        x[i,0]=nu;x[j,0]=nv
    return x


@njit(cache=True)
def transform(matrix,v):
    out=np.zeros((3,v.shape[1],2))
    for i in range(3):
        for j in range(v.shape[1]):
            q=(0.,0.)
            for k in range(3):
                q=plus(q,times(matrix[i,k],v[k,j]))
            out[i,j]=q
    return out


@njit(cache=True)
def intersect_mean(jet,center,rad):
    for i in range(3):
        error=0.
        for j in range(3):
            v=maximum_abs(jet[i,j+1])
            inc=times((v,v),(rad[j],rad[j]))[1]
            error=np.nextafter(error+inc,np.inf)
        jet[i,0]=overlap(jet[i,0],plus(center[i,0],(-error,error)))


@njit(cache=True)
def leaf(mid,rad,order,p,base,origin,inv,pi2,cs,eps):
    center=np.zeros((3,1,2));wide=np.zeros((3,1,2))
    for i in range(3):
        center[i,0]=(mid[i],mid[i]);wide[i,0]=plus(center[i,0],(-rad[i],rad[i]))
    center=transform(base,center);wide=transform(base,wide)
    jet=np.empty((3,4,2))
    for i in range(3):
        center[i,0]=plus(center[i,0],origin[i]);jet[i,0]=plus(wide[i,0],origin[i])
        for j in range(3):
            jet[i,j+1]=base[i,j]
    for k in range(1,order+1):
        center=rotate(center,p,pi2,cs,eps);jet=rotate(jet,p,pi2,cs,eps)
        intersect_mean(jet,center,rad)
        ci=center.copy();ji=jet.copy()
        for i in range(3):
            ci[i,0]=minus(ci[i,0],origin[i]);ji[i,0]=minus(ji[i,0],origin[i])
        ci=transform(inv,ci);ji=transform(inv,ji)
        intersect_mean(ji,ci,rad)
        if k<order:
            if not any_outside(ji):
                return False
        else:
            for i in range(3):
                if not (ji[i,0,0]>-1 and ji[i,0,1]<1):
                    return False
    return True


@njit(inline="always")
def any_outside(jet):
    for i in range(3):
        if jet[i,0,0]>1 or jet[i,0,1]<-1:
            return True
    return False


def frame(packet):
    q=Fraction.from_float
    b=[[q(v) for v in row] for row in packet["basis"][0]]
    l=list(map(q,packet["lo"][0]));h=list(map(q,packet["hi"][0]))
    origin=[q(x)+sum(b[i][j]*(h[j]+l[j])/2 for j in range(3)) for i,x in enumerate(packet["center"][0])]
    mat=[[b[i][j]*(h[j]-l[j])/2 for j in range(3)] for i in range(3)]
    def minor(i,j):
        a=[a for a in range(3) if a!=i];c=[a for a in range(3) if a!=j]
        return mat[a[0]][c[0]]*mat[a[1]][c[1]]-mat[a[0]][c[1]]*mat[a[1]][c[0]]
    det=sum((-1)**j*mat[0][j]*minor(0,j) for j in range(3))
    inverse=[[(-1)**(i+j)*minor(j,i)/det for j in range(3)] for i in range(3)]
    return (np.array([[enclosure(x) for x in row] for row in mat]),np.array([enclosure(x) for x in origin]),
            np.array([[enclosure(x) for x in row] for row in inverse]))


def audit(box,proof,parameter_bits=None):
    order=box["case_order"]
    if parameter_bits is None:
        p=np.array([[v,v] for v in box["params"]])
    else:
        radius=Fraction(1,2**parameter_bits)
        p=np.array([[enclosure(Fraction.from_float(v)-radius)[0],enclosure(Fraction.from_float(v)+radius)[1]]
                    for v in box["params"]])
    if not (p[0,0]>0 and p[0,1]<1):
        raise ValueError("Need a strictly positive contraction interval below one")
    m,o,iv=frame(box);pi2,cs,eps=trig_constants()
    tree=zlib.decompress(base64.b64decode(proof["partition_zlib_base64"]))
    entire=(0,0,0);stack=[(entire,entire,True)]
    pos=accepted=extra=0
    while stack:
        indices,depths,frozen=stack.pop()
        if max(depths)>48:
            raise RuntimeError("Beyond exact binary64 dyadic depth")
        axis=3
        if frozen:
            axis=tree[pos];pos+=1
            if axis>3:
                raise ValueError("Malformed partition")
        if axis==3:
            rad=np.array([2.**-d for d in depths])
            mid=np.array([(2*i+1)*h-1 for i,h in zip(indices,rad)])
            if leaf(mid,rad,order,p,m,o,iv,pi2,cs,eps):
                accepted+=1;continue
            axis=int(np.argmax(rad));frozen=False;extra+=1
            if extra>1_000_000:
                raise RuntimeError("Extra subdivision budget exhausted")
        for sign in (0,1):
            ii=list(indices);dd=list(depths)
            ii[axis]=2*ii[axis]+sign;dd[axis]+=1
            stack.append((tuple(ii),tuple(dd),frozen))
    if pos!=len(tree):
        raise RuntimeError("Trailing partition data")
    return dict(case=box["case"],phase_regions=order,stored_tree_nodes_checked=pos,accepted_leaves=accepted,
                extra_splits=extra,independent_parameter_radius_bits=parameter_bits,
                method="IEEE binary64 outward intervals, rational frame and Taylor trigonometry",
                pi_over_two_interval=pi2,reduced_angle_remainder_upper=eps)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("box",type=Path);ap.add_argument("certificate",type=Path)
    ap.add_argument("--output",type=Path)
    ap.add_argument("--parameter-radius-bits",type=int)
    a=ap.parse_args()
    result=audit(json.loads(a.box.read_text()),json.loads(a.certificate.read_text()),a.parameter_radius_bits)
    if a.output:
        a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)

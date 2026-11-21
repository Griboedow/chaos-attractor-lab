"""Independent mpmath interval replay, exact rational frames and forward-mode AD."""
import argparse
import base64
import json
import zlib
from fractions import Fraction
from pathlib import Path
from mpmath import iv


def inverse_rational(a):
    rows=[a[i].copy()+[Fraction(i==j) for j in range(3)] for i in range(3)]
    for j in range(3):
        i=next(i for i in range(j,3) if rows[i][j])
        rows[j],rows[i]=rows[i],rows[j]
        x=rows[j][j]
        rows[j]=[v/x for v in rows[j]]
        for i in range(3):
            if i!=j:
                x=rows[i][j]
                rows[i]=[u-x*v for u,v in zip(rows[i],rows[j])]
    return [r[3:] for r in rows]


def exact_frame(packet):
    q=Fraction.from_float
    c=list(map(q,packet["center"][0]))
    b=[list(map(q,r)) for r in packet["basis"][0]]
    l=list(map(q,packet["lo"][0]));h=list(map(q,packet["hi"][0]))
    width=[(y-x)/2 for x,y in zip(l,h)]
    mid=[(y+x)/2 for x,y in zip(l,h)]
    o=[c[i]+sum(row[j]*mid[j] for j in range(3)) for i,row in enumerate(b)]
    m=[[v*width[j] for j,v in enumerate(r)] for r in b]
    convert=lambda x:iv.mpf(x.numerator)/x.denominator
    return list(map(convert,o)),[[convert(v) for v in r] for r in m],[[convert(v) for v in r] for r in inverse_rational(m)]


class Jet:
    __slots__=("x","d")
    def __init__(self,x,d):
        self.x=x;self.d=d

    def __add__(self,other):
        if isinstance(other,Jet):
            return Jet(self.x+other.x,[a+b for a,b in zip(self.d,other.d)])
        return Jet(self.x+other,self.d)
    __radd__=__add__

    def __sub__(self,other):
        return Jet(self.x-other.x,[a-b for a,b in zip(self.d,other.d)])

    def __mul__(self,other):
        if isinstance(other,Jet):
            return Jet(self.x*other.x,[a*other.x+self.x*b for a,b in zip(self.d,other.d)])
        return Jet(self.x*other,[a*other for a in self.d])
    __rmul__=__mul__

    def __pow__(self,n):
        return Jet(self.x**n,[n*self.x**(n-1)*a for a in self.d])


def sincos(t):
    if not isinstance(t,Jet):
        return iv.sin(t),iv.cos(t)
    si,co=iv.sin(t.x),iv.cos(t.x)
    return Jet(si,[co*v for v in t.d]),Jet(co,[-si*v for v in t.d])


def recurrence(state,p):
    x,y,z=state
    x,y,z=x*p[0]+(1-p[0]),y*p[0],z*p[0]
    x,y=pair(x,y,z*p[4]+z**3*p[7]+p[1])
    y,z=pair(y,z,x*p[5]+x**3*p[8]+p[2])
    z,x=pair(z,x,y*p[6]+y**3*p[9]+p[3])
    return x,y,z


def pair(x,y,angle):
    si,co=sincos(angle)
    u,v=x*co-y*si,x*si+y*co
    if isinstance(angle,Jet):
        # The product-rule evaluation and the rotating-frame identity enclose
        # the same derivative; intersection controls dependency overestimation.
        u.d=[intersect(w,co.x*a-si.x*b-v.x*c) for w,a,b,c in zip(u.d,x.d,y.d,angle.d)]
        v.d=[intersect(w,si.x*a+co.x*b+u.x*c) for w,a,b,c in zip(v.d,x.d,y.d,angle.d)]
    return u,v


def transform(m,x):
    return [x[0]*r[0]+x[1]*r[1]+x[2]*r[2] for r in m]


def intersect(a,b):
    lo=max(a.a,b.a);hi=min(a.b,b.b)
    if not lo<=hi:
        raise RuntimeError("Disjoint enclosures")
    return iv.mpf([lo,hi])


def limit(center,jet,radius):
    error=sum(abs(v)*d for v,d in zip(jet.d,radius))
    return intersect(jet.x,center+iv.mpf([-error.b,error.b]))


def splits(mid,rad,axis):
    half=rad.copy();half[axis]/=2
    for sign in (-1,1):
        c=mid.copy();c[axis]+=sign*half[axis]
        if Fraction.from_float(c[axis])-sign*Fraction.from_float(half[axis])!=Fraction.from_float(mid[axis]):
            raise RuntimeError("Inexact dyadic child")
        yield c,half


def replay(packet,cert,digits,limit_extra):
    iv.dps=digits
    org,base,inv=exact_frame(packet)
    params=[iv.mpf(v) for v in packet["params"]]
    tree=zlib.decompress(base64.b64decode(cert["partition_zlib_base64"]))
    pos=0;accepted=extra=0
    stack=[([0.,0.,0.],[1.,1.,1.],True)]
    while stack:
        mid,rad,original=stack.pop()
        if original:
            if pos>=len(tree):
                raise RuntimeError("Truncated subdivision tree")
            axis=tree[pos];pos+=1
            if axis<3:
                stack.extend((c,h,True) for c,h in splits(mid,rad,axis))
                continue
            if axis!=3:
                raise ValueError("Invalid tree entry")
        spans=[iv.mpf([m-r,m+r]) for m,r in zip(mid,rad)]
        center=[a+b for a,b in zip(org,transform(base,mid))]
        state=[Jet(a+b,base[i]) for i,(a,b) in enumerate(zip(org,transform(base,spans)))]
        passed=True
        for k in range(1,packet["case_order"]+1):
            center=recurrence(center,params);state=recurrence(state,params)
            state=[Jet(limit(a,b,rad),b.d) for a,b in zip(center,state)]
            mc=transform(inv,[a-b for a,b in zip(center,org)])
            jm=transform(inv,[a+(-b) for a,b in zip(state,org)])
            bound=[limit(a,b,rad) for a,b in zip(mc,jm)]
            if k<packet["case_order"]:
                passed=any(v.a>1 or v.b< -1 for v in bound)
            else:
                passed=all(v.a> -1 and v.b<1 for v in bound)
            if not passed:
                break
        if passed:
            accepted+=1
        else:
            extra+=1
            if extra>limit_extra:
                raise RuntimeError("Exceeded secondary subdivision allowance")
            allowed=[j for j in range(3) if rad[j]>=max(rad)/64]
            axis=max(allowed,key=lambda j:rad[j]*max(float(abs(b.d[j]).b) for b in jm))
            stack.extend((c,h,False) for c,h in splits(mid,rad,axis))
    if pos!=len(tree):
        raise RuntimeError("Trailing subdivision nodes")
    entry=None
    if cert.get("primary_start_entry"):
        record=cert["primary_start_entry"]
        iv.dps=max(digits,int(record["precision_bits"]*.302)+20)
        org,_,inv=exact_frame(packet)
        x=[iv.mpf(v) for v in (.1,.2,.3)]
        for _ in range(record["iterate"]):
            x=recurrence(x,params)
        coords=transform(inv,[a-b for a,b in zip(x,org)])
        if not all(v.a> -1 and v.b<1 for v in coords):
            raise RuntimeError("Independent start-entry replay failed")
        entry=dict(iterate=record["iterate"],decimal_digits=iv.dps)
    return dict(case=packet["case"],phase_regions=packet["case_order"],interval_decimal_digits=digits,
                stored_tree_nodes_checked=pos,accepted_leaf_boxes=accepted,extra_splits=extra,
                independently_verified_start_entry=entry,
                method="mpmath.iv with forward-mode automatic differentiation and rational inverse frames")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("box",type=Path);ap.add_argument("certificate",type=Path)
    ap.add_argument("--digits",type=int,default=45)
    ap.add_argument("--extra",type=int,default=10000)
    ap.add_argument("--output",type=Path)
    a=ap.parse_args()
    r=replay(json.loads(a.box.read_text()),json.loads(a.certificate.read_text()),a.digits,a.extra)
    if a.output:
        a.output.write_text(json.dumps(r,indent=2)+"\n")
    print(r,flush=True)

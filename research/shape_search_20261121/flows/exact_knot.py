"""Exact crossing combinatorics of an integer polygon in an integer projection."""
from bisect import bisect_left
from fractions import Fraction
import numpy as np
from knot_diagram import diagram_invariants


def exact_diagram(vertices,matrix,gauss=False):
    # Convert to Python integers before all orientation and depth predicates.
    data=[[sum(int(p[k])*int(matrix[k][j]) for k in range(3)) for j in range(3)] for p in vertices]
    n=len(data)
    edges=[]
    for i,p in enumerate(data):
        q=data[(i+1)%n]
        edges.append((p,tuple(b-a for a,b in zip(p,q)),min(p[0],q[0]),max(p[0],q[0]),min(p[1],q[1]),max(p[1],q[1])))
    crossings=[]
    for i,(p,a,x0,x1,y0,y1) in enumerate(edges):
        for j in range(i+2,n):
            if i==0 and j==n-1:
                continue
            q,b,u0,u1,v0,v1=edges[j]
            if x1<u0 or u1<x0 or y1<v0 or v1<y0:
                continue
            dx=q[0]-p[0];dy=q[1]-p[1]
            den=a[0]*b[1]-a[1]*b[0]
            sn=dx*b[1]-dy*b[0];tn=dx*a[1]-dy*a[0]
            if den==0:
                if sn==tn==0:
                    raise ValueError('Collinear projected segments')
                continue
            orient=1 if den>0 else -1
            if den<0:
                den=-den;sn=-sn;tn=-tn
            if sn<0 or sn>den or tn<0 or tn>den:
                continue
            if sn in (0,den) or tn in (0,den):
                raise ValueError('Projected crossing at polygon vertex')
            depth=(p[2]-q[2])*den+sn*a[2]-tn*b[2]
            if not depth:
                raise ValueError('Intersecting spatial edges')
            at_i=Fraction(i*den+sn,den)
            at_j=Fraction(j*den+tn,den)
            crossings.append((at_j,at_i,orient) if depth>0 else (at_i,at_j,-orient))
    if gauss:
        events=[]
        for idx,(under,over,sign) in enumerate(crossings):
            events.extend([(under,idx,False,sign),(over,idx,True,sign)])
        events.sort()
        if len(set(v[0] for v in events))!=len(events):raise ValueError('Repeated crossing coordinate')
        return [(idx,over,sign) for _,idx,over,sign in events]
    if not crossings:
        return dict(crossings=0,determinant=1,coefficients=[1])
    crossings.sort()
    all_positions=[v for u,o,_ in crossings for v in (u,o)]
    if len(set(all_positions))!=len(all_positions):
        raise ValueError('Repeated crossing coordinate')
    sites=[u for u,_,_ in crossings]
    count=len(sites)
    arcs=[((bisect_left(sites,o)-1)%count,(i-1)%count,i,sign) for i,(_,o,sign) in enumerate(crossings)]
    return diagram_invariants(arcs,None,True)


def prepare(cloud,camera):
    shift=cloud.mean(0)
    scale=np.sqrt(np.mean(np.sum((cloud-shift)**2,axis=1)))
    integer=np.rint((cloud-shift)/scale*10**9).astype(np.int64)
    rounded=np.rint(np.linalg.qr(np.random.default_rng(camera).normal(size=(3,3)))[0]*10000).astype(int)
    # An arbitrary invertible integer coordinate change preserves the polynomial.
    a,b,c=rounded.tolist()
    det=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
    if det==0:
        raise ValueError('Singular camera')
    return integer,rounded,shift,scale

"""Integer Alexander invariants of generic polygon projections."""
from fractions import Fraction
import numpy as np
from numba import njit
try:
    from flint import fmpz_mat
except ImportError:
    fmpz_mat=None


@njit(cache=True)
def projected_crossings(z):
    n=len(z)
    under=[]
    over=[]
    orientation=[]
    minimum_depth=1e30
    for i in range(n):
        i1=(i+1)%n
        px,py,pz=z[i]; ax,ay,az=z[i1]-z[i]
        for j in range(i+2,n):
            if i==0 and j==n-1:
                continue
            j1=(j+1)%n
            qx,qy,qz=z[j]; bx,by,bz=z[j1]-z[j]
            den=ax*by-ay*bx
            if abs(den)<1e-15:
                continue
            s=((qx-px)*by-(qy-py)*bx)/den
            t=((qx-px)*ay-(qy-py)*ax)/den
            if 0<s<1 and 0<t<1:
                dz=(pz+s*az)-(qz+t*bz)
                minimum_depth=min(minimum_depth,abs(dz))
                if dz>0:
                    over.append(i+s);under.append(j+t);orientation.append(1 if den>0 else -1)
                else:
                    under.append(i+s);over.append(j+t);orientation.append(-1 if den>0 else 1)
    return np.array(under),np.array(over),np.array(orientation),minimum_depth


def determinant(a):
    """Fraction-free Gaussian elimination on an integer matrix."""
    size=len(a)
    if size==0:
        return 1
    if size>12 and fmpz_mat is not None:
        return int(fmpz_mat(a).det())
    sign=1;previous=1
    for k in range(size-1):
        if not a[k][k]:
            pivot=next((i for i in range(k+1,size) if a[i][k]),None)
            if pivot is None:
                return 0
            a[k],a[pivot]=a[pivot],a[k];sign*=-1
        diagonal=a[k][k]
        for i in range(k+1,size):
            for j in range(k+1,size):
                a[i][j]=(a[i][j]*diagonal-a[i][k]*a[k][j])//previous
        previous=diagonal
    return sign*a[-1][-1]


def relation_matrix(arcs, t):
    n=len(arcs)
    mat=[[0]*n for _ in range(n)]
    for i,(over,entry,exit,sign) in enumerate(arcs):
        mat[i][over]+=1-t
        if sign<0:
            entry,exit=exit,entry
        mat[i][entry]+=t
        mat[i][exit]-=1
    return [row[:-1] for row in mat[:-1]]


def alexander(z, seed, polynomial=False):
    random=np.random.default_rng(seed)
    rotation=np.linalg.qr(random.normal(size=(3,3)))[0]
    view=z@rotation
    under,over,orientation,depth=projected_crossings(view)
    n=len(under)
    if n>100:
        return None
    if n==0:
        return dict(crossings=0,determinant=1,coefficients=[1],minimum_depth=None)
    order=np.argsort(under)
    sites=under[order]
    arcs=[((np.searchsorted(sites,over[i])-1)%n,(j-1)%n,j,int(orientation[i])) for j,i in enumerate(order)]
    return diagram_invariants(arcs,float(depth),polynomial)


def diagram_invariants(arcs,depth,polynomial):
    n=len(arcs)
    det=abs(determinant(relation_matrix(arcs,-1)))
    report=dict(crossings=n,determinant=det,minimum_depth=depth)
    if not polynomial:
        return report
    values=[determinant(relation_matrix(arcs,t)) for t in range(n)]
    basis=[Fraction(1)]
    coefficients=[Fraction(0)]*n
    for k in range(n):
        for j,v in enumerate(basis):
            coefficients[j]+=v*values[0]
        values=[right-left for left,right in zip(values[:-1],values[1:])]
        nxt=[Fraction(0)]*(len(basis)+1)
        for j,v in enumerate(basis):
            nxt[j]-=k*v/(k+1);nxt[j+1]+=v/(k+1)
        basis=nxt
    while len(coefficients)>1 and coefficients[-1]==0:
        coefficients.pop()
    while len(coefficients)>1 and coefficients[0]==0:
        coefficients.pop(0)
    if any(v.denominator!=1 for v in coefficients):
        raise RuntimeError('Nonintegral Alexander coefficients')
    sign=1 if coefficients[0]>0 else -1
    a=[sign*v.numerator for v in coefficients]
    if a!=a[::-1] or abs(sum(a))!=1 or abs(sum((-1)**j*v for j,v in enumerate(a)))!=det:
        raise RuntimeError('Alexander normalization did not pass')
    report['coefficients']=a
    return report


def positive_controls():
    t=np.linspace(0,2*np.pi,500,endpoint=False)+.00013
    ring=np.column_stack((np.cos(t),np.sin(t),np.zeros_like(t)))
    def torus(p,q):
        return np.column_stack(((2+np.cos(q*t))*np.cos(p*t),(2+np.cos(q*t))*np.sin(p*t),np.sin(q*t)))
    eight=np.column_stack(((2+np.cos(2*t))*np.cos(3*t),(2+np.cos(2*t))*np.sin(3*t),np.sin(4*t)))
    for z,expected in ((ring,[1]),(torus(2,3),[1,-1,1]),(torus(2,5),[1,-1,1,-1,1]),(eight,[1,-3,1])):
        for seed in [17,35,81]:
            assert alexander(z,seed,True)['coefficients']==expected


if __name__=='__main__':
    positive_controls()
    print('Circle, trefoil, cinquefoil and figure-eight controls passed in three projections.')

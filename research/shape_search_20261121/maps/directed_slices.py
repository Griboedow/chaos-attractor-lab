"""Independently bound inverse distances between phase points with native intervals."""
import argparse
import json
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from directed_return import enclosure,frame,intersect_mean,maximum_abs,minus,plus,rotate,sincos,times,transform,trig_constants
import numpy as np
from numba import njit


@njit(cache=True)
def phase_hulls(M,o,p,steps,grid,pi2,cs,eps):
    n=2**grid
    h=1./n
    rad=np.array([h,h,h])
    hull=np.empty((steps,3,2))
    hull[:,:,0]=np.inf
    hull[:,:,1]=-np.inf
    for a in range(n):
        for b in range(n):
            for c in range(n):
                mid=np.array([(2*a+1)*h-1,(2*b+1)*h-1,(2*c+1)*h-1])
                center=np.zeros((3,1,2));wide=center.copy()
                for i in range(3):
                    center[i,0]=(mid[i],mid[i]);wide[i,0]=plus(center[i,0],(-h,h))
                center=transform(M,center);wide=transform(M,wide)
                jet=np.zeros((3,4,2))
                for i in range(3):
                    center[i,0]=plus(center[i,0],o[i]);jet[i,0]=plus(wide[i,0],o[i])
                    for j in range(3):
                        jet[i,j+1]=M[i,j]
                for j in range(steps):
                    for k in range(3):
                        hull[j,k,0]=min(hull[j,k,0],jet[k,0,0])
                        hull[j,k,1]=max(hull[j,k,1],jet[k,0,1])
                    center=rotate(center,p,pi2,cs,eps,True)
                    jet=rotate(jet,p,pi2,cs,eps,True)
                    intersect_mean(jet,center,rad)
    return hull


@njit(cache=True)
def inverse(x,p,pi2,cs,eps):
    J=np.zeros((3,3,2))
    for i in range(3):
        J[i,i]=(1.,1.)
    for i in (2,1,0):
        j,k=(i+1)%3,(i+2)%3
        w=x[k].copy();u=x[i].copy();v=x[j].copy()
        sq=times(w,w)
        angle=minus((0.,0.),plus(p[i+1],times(w,plus(p[i+4],times(sq,p[i+7])))))
        slope=minus((0.,0.),plus(p[i+4],times((3.,3.),times(sq,p[i+7]))))
        si,co=sincos(angle,pi2,cs,eps)
        x[i]=minus(times(co,u),times(si,v));x[j]=plus(times(si,u),times(co,v))
        for l in range(3):
            du,dv=J[i,l].copy(),J[j,l].copy()
            dangle=times(slope,J[k,l])
            J[i,l]=minus(minus(times(co,du),times(si,dv)),times(x[j],dangle))
            J[j,l]=plus(plus(times(si,du),times(co,dv)),times(x[i],dangle))
    reciprocal=(np.nextafter(1./p[0,1],-np.inf),np.nextafter(1./p[0,0],np.inf))
    for i in range(3):
        for j in range(3):
            J[i,j]=times(J[i,j],reciprocal)
    return J


@njit(cache=True)
def derivative_hull(native,p,grid,pi2,cs,eps):
    n=2**grid;h=1./n
    hull=np.empty((3,3,2));hull[:,:,0]=np.inf;hull[:,:,1]=-np.inf
    for a in range(n):
        for b in range(n):
            for c in range(n):
                index=(a,b,c)
                x=np.empty((3,2))
                for i in range(3):
                    lo,hi=native[i,0],native[i,1]
                    diff=minus((hi,hi),(lo,lo))
                    left=index[i]*h;right=(index[i]+1)*h
                    x[i]=plus((lo,lo),times(diff,(left,right)))
                mat=inverse(x,p,pi2,cs,eps)
                for i in range(3):
                    for j in range(3):
                        hull[i,j,0]=min(hull[i,j,0],mat[i,j,0])
                        hull[i,j,1]=max(hull[i,j,1],mat[i,j,1])
    return hull


def verify(box,proof,grid):
    m=proof["phase_regions"];bits=proof["parameter_radius_bits"]
    if m!=box["case_order"] or proof["case"]!=box["case"]:
        raise ValueError("Wrong cyclic proof")
    epsp=Fraction(1,2**bits)
    p=np.array([[enclosure(Fraction.from_float(v)-epsp)[0],enclosure(Fraction.from_float(v)+epsp)[1]] for v in box["params"]])
    M,o,_=frame(box);pi2,cs,eps=trig_constants()
    hull=phase_hulls(M,o,p,m,grid,pi2,cs,eps)
    matrices=[derivative_hull(v,p,grid,pi2,cs,eps) for v in hull]
    J=np.zeros((3,3,2))
    for i in range(3):
        J[i,i]=[1.,1.]
    for j in (0,*reversed(range(1,m))):
        J=transform(matrices[j],J)
    ctx.prec=128
    L=sum(arb(float(maximum_abs(J[i,j])))**2 for i in range(3) for j in range(3)).sqrt()
    low=1+arb(2).log()/L.log()
    q=Fraction(int(float(low.lower())*1e6),1_000_000)
    if not low>arb(q.numerator)/q.denominator:
        raise ArithmeticError("Rounding failure")
    return dict(case=box["case"],parameter_radius_bits=bits,grid_cells_per_region=2**(3*grid),
        independent_method="outward binary64 enclosure with rational trigonometry",
        native_phase_hulls=hull.tolist(),inverse_return_Lipschitz_upper_approximation=float(L.upper()),
        rigorous_dimension_lower_rational=str(q),dimension_lower=float(q))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("box",type=Path);ap.add_argument("trap",type=Path)
    ap.add_argument("--grid-power",type=int,default=5);ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    result=verify(json.loads(a.box.read_text()),json.loads(a.trap.read_text()),a.grid_power)
    a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(result["case"],result["dimension_lower"],flush=True)

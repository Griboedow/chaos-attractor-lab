"""Tighten the inverse Lipschitz estimate on the five certified phase neighborhoods."""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,arb_mat,ctx
from slice_bound import ball
from verify_boxes import column,frames,rotate,values
from verify_return import tighten


def hull(v,w):
    return w if v is None else v.union(w)


def endpoints(x):
    lo=math.nextafter(float(x.lower()),-math.inf)
    hi=math.nextafter(float(x.upper()),math.inf)
    if not (arb(lo)<=x and arb(hi)>=x):
        raise ArithmeticError("Endpoint rounding failed")
    return (lo,hi)


def affine_frame(bounds):
    mid=[(arb(lo)+arb(hi))/2 for lo,hi in bounds]
    rad=[(arb(hi)-arb(lo))/2 for lo,hi in bounds]
    return mid,rad


def cubes(power):
    n=2**power
    width=1/n
    for index in itertools.product(range(n),repeat=3):
        yield [(2*j+1)*width-1 for j in index],[width]*3


def inverse_jac(x,p,identity):
    J=identity
    x=list(x)
    for i in (2,1,0):
        j,k=(i+1)%3,(i+2)%3
        w=x[k]
        angle=-(p[1+i]+p[4+i]*w+p[7+i]*w**3)
        slope=-(p[4+i]+3*p[7+i]*w*w)
        s,c=angle.sin_cos()
        u,v=c*x[i]-s*x[j],s*x[i]+c*x[j]
        rot=arb_mat([[int(a==b) for b in range(3)] for a in range(3)])
        rot[i,i],rot[i,j],rot[i,k]=c,-s,-v*slope
        rot[j,i],rot[j,j],rot[j,k]=s,c,u*slope
        J=rot*J
        x[i],x[j]=u,v
    return J*(1/p[0])


def norm(A):
    return sum(A[i,j].abs_upper()**2 for i in range(3) for j in range(3)).sqrt()


def compute(box,symbolic,trap,grid):
    if symbolic["parameter_radius_bits"]!=trap["parameter_radius_bits"] or not symbolic["all_faces_and_full_slabs_verified"]:
        raise ValueError("Different parameter blocks or missing symbolic proof")
    if any(a["case"]!=box["case"] for a in (symbolic,trap)):
        raise ValueError("Different maps")
    ctx.prec=128
    bits=symbolic["parameter_radius_bits"]
    p=[arb(float(x),2.**-bits) for x in box["params"]]
    m=trap["phase_regions"]
    if m!=symbolic["cover"]["map_steps"] or m!=box["case_order"]:
        raise ValueError("Different return maps")
    o,M,_=frames(box)
    o,M=o[0],M[0]
    boxes=[[None]*3 for _ in range(m)]
    for mid,rad in cubes(grid):
        center=o+M*column(mid)
        wide=o+M*column([arb(v,r) for v,r in zip(mid,rad)])
        J=M
        for i in range(m):
            for j in range(3):
                boxes[i][j]=hull(boxes[i][j],wide[j,0])
            if i+1<m:
                center,_=rotate(values(center),p)
                wide,J=rotate(values(wide),p,J)
                wide=tighten(wide,center,J,rad)
    identity=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    bounds=[[endpoints(x) for x in row] for row in boxes]
    row_norms=[]
    derivative=[]
    for box_bounds in bounds:
        center,rad=affine_frame(box_bounds)
        upper=arb(0)
        enclosure=[[None]*3 for _ in range(3)]
        for local_mid,local_rad in cubes(grid):
            x=[v+h*arb(w,r) for v,h,w,r in zip(center,rad,local_mid,local_rad)]
            J=inverse_jac(x,p,identity)
            upper=max(upper,norm(J).abs_upper())
            for j in range(3):
                for k in range(3):
                    enclosure[j][k]=hull(enclosure[j][k],J[j,k])
        row_norms.append(upper)
        derivative.append(arb_mat(enclosure))
    product=identity
    for j in (0,*reversed(range(1,m))):
        product=derivative[j]*product
    product_bound=math.prod(row_norms)
    matrix_bound=norm(product)
    L=min(product_bound.abs_upper(),matrix_bound.abs_upper())
    if not L>1:
        raise ArithmeticError("Invalid inverse expansion bound")
    low=1+arb(2).log()/L.log()
    rounded=Fraction(int(float(low.lower())*1_000_000),1_000_000)
    if not low>ball(rounded):
        raise ArithmeticError("Decimal lower bound did not round down")
    return dict(case=box["case"],parameter_radius_bits=bits,grid_cells_per_region=2**(3*grid),
        precision_bits=ctx.prec,phase_enclosures_binary64=bounds,
        phase_inverse_Lipschitz_upper_approximations=[float(v) for v in row_norms],
        inverse_return_product_of_norms_upper_approximation=float(product_bound.upper()),
        inverse_return_interval_matrix_norm_upper_approximation=float(matrix_bound.upper()),
        rigorous_dimension_lower_rational=str(rounded),dimension_lower=float(rounded))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for arg in ("box","symbolic","trap"):
        ap.add_argument(arg,type=Path)
    ap.add_argument("--grid-power",type=int,default=3)
    ap.add_argument("--output",required=True,type=Path)
    a=ap.parse_args()
    result=compute(*(json.loads(getattr(a,k).read_text()) for k in ("box","symbolic","trap")),a.grid_power)
    a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(result["case"],result["dimension_lower"],flush=True)

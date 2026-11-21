"""Bound Hausdorff dimension on a certified return box by singular-value contraction."""
import argparse
import base64
import hashlib
import json
import zlib
from pathlib import Path
from flint import arb, arb_mat, ctx
from verify_boxes import column, frames, rotate, values
from verify_return import tighten


def exterior_step(x,p,area):
    # In 3-space the action on oriented area vectors is det(J) J^(-T).
    # For a twist in i,j, the k-row of its area action is
    # [old_y * angle', -old_x * angle', 1] in the ordered coordinates i,j,k.
    x=[p[0]*v for v in x]
    x[0]+=1-p[0]
    area=area*(p[0]*p[0])
    for i in range(3):
        j,k=(i+1)%3,(i+2)%3
        u,v,w=x[i],x[j],x[k]
        slope=p[4+i]+3*p[7+i]*w*w
        s,c=(p[1+i]+p[4+i]*w+p[7+i]*w**3).sin_cos()
        block=arb_mat([[int(a==b) for b in range(3)] for a in range(3)])
        block[i,i],block[i,j],block[j,i],block[j,j]=c,-s,s,c
        block[k,i],block[k,j]=v*slope,-u*slope
        area=block*area
        x[i],x[j]=c*u-s*v,s*u+c*v
    return area


def cofactors(J):
    return arb_mat([[J[(i+1)%3,(j+1)%3]*J[(i+2)%3,(j+2)%3]
                    -J[(i+1)%3,(j+2)%3]*J[(i+2)%3,(j+1)%3] for j in range(3)] for i in range(3)])


def frobenius_squared(a):
    return sum(a[i,j].abs_upper()**2 for i in range(3) for j in range(3))


def verify(box,num,den,budget,bits):
    if not (0<num<den):
        raise ValueError("Use a dimension strictly between one and two")
    ctx.prec=128
    o,M,_=frames(box);o,M=o[0],M[0]
    ident=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    p=[arb(float(v)) if bits is None else arb(float(v),2.**-bits) for v in box["params"]]
    stack=[([0.]*3,[1.]*3)]
    nodes=leaves=0
    accepted_max=arb(0)
    encoding=bytearray()
    while stack:
        mid,rad=stack.pop()
        nodes+=1
        if nodes>budget:
            raise RuntimeError("Dimension subdivision budget exhausted")
        center=o+M*column(mid)
        wide=o+M*column([arb(c,r) for c,r in zip(mid,rad)])
        J=ident
        C=ident
        for _ in range(box["case_order"]):
            C=exterior_step(values(wide),p,C)
            center,_=rotate(values(center),p)
            wide,J=rotate(values(wide),p,J)
            wide=tighten(wide,center,J*M,rad)
        # Compute cofactor bounds in two ways and intersect each entry.
        cross=cofactors(J)
        C=arb_mat([[C[i,j].intersection(cross[i,j]) for j in range(3)] for i in range(3)])
        upper=frobenius_squared(J)**(den-num)*frobenius_squared(C)**num
        if upper<1:
            encoding.append(3)
            accepted_max=max(accepted_max,upper.abs_upper());leaves+=1
            continue
        physical=[sum(abs(M[i,j]).abs_upper() for i in range(3))*rad[j] for j in range(3)]
        axis=max(range(3),key=lambda j:float(physical[j]))
        dr=rad.copy();dr[axis]/=2
        if not dr[axis] or dr[axis]<2**-43:
            raise RuntimeError("Unresolved point at numerical partition depth")
        encoding.append(axis)
        for sign in (-1,1):
            m=mid.copy();m[axis]+=sign*dr[axis]
            if arb(m[axis])-sign*arb(dr[axis])!=arb(mid[axis]):
                raise ValueError("Inexact split")
            stack.append((m,dr))
    return dict(case=box["case"],map_steps=box["case_order"],parameter_radius_bits=bits,
        dimension_upper_rational=f"{den+num}/{den}",dimension_upper=(den+num)/den,
        power_used=2*den,maximum_accepted_power_approximation=float(accepted_max),
        nodes=nodes,leaves=leaves,precision_bits=ctx.prec,
        box_digest=hashlib.sha256(json.dumps(box,sort_keys=True).encode()).hexdigest(),
        partition_zlib_base64=base64.b64encode(zlib.compress(encoding,9)).decode("ascii"))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("box",type=Path)
    ap.add_argument("--fraction",type=int,nargs=2,required=True,metavar=("NUM","DEN"))
    ap.add_argument("--budget",type=int,default=200_000)
    ap.add_argument("--parameter-radius-bits",type=int)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    result=verify(json.loads(a.box.read_text()),*a.fraction,a.budget,a.parameter_radius_bits)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(result["case"],result["dimension_upper_rational"],result["nodes"],flush=True)

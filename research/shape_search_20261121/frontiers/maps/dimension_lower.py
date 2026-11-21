"""Evaluate exact shear bounds and the slice gap throughout a coefficient block."""
import argparse
import json
from fractions import Fraction
from pathlib import Path
from entropy import FACES,INPUT
from flint import arb,ctx
from verify_boxes import column,frames,values
from slice_bound import ball

def bound(box,bits):
    ctx.prec=192
    o,M,V=frames(box)
    cuboid=o[0]+M[0]*column([arb(0,1)]*3)
    radius_squared=sum(v.abs_upper()**2 for v in values(cuboid))
    if not radius_squared<1:raise ArithmeticError("The source box must be in the unit ball")
    pars=[Fraction.from_float(x) for x in box["params"]]
    eps=Fraction(1,2**bits) if bits else Fraction(0)
    if not 0<pars[0]-eps<pars[0]+eps<1 and not (eps==0 and 0<pars[0]<1):
        raise ValueError("Affine contraction range")
    L=1/ball(pars[0]-eps);terms=[]
    for i in range(3):
        k=abs(pars[4+i])+eps;h=3*(abs(pars[7+i])+eps)
        ys=[Fraction(0),Fraction(1)]
        if h:
            v=(2*h-k)/(3*h)
            if 0<v<1:ys.append(v)
        q=max((k+h*y)**2*(1-y) for y in ys)
        a=ball(q).sqrt();L*=((4+a**2).sqrt()+a)/2
        terms.append(str(q))
    N=2*box["case_order"]
    invG=L**N;d=1+arb(2).log()/invG.log()
    lower=Fraction(int(float(d.lower())*10**6),10**6)
    grad=sum(V[0][2,j]**2 for j in range(3)).sqrt()
    gap=(arb(FACES[2])-arb(FACES[1]))/grad
    gap_lower=Fraction(int(float(gap.lower())*10**9),10**9)
    if not (d>ball(lower) and gap>ball(gap_lower) and gap_lower>0):raise ArithmeticError("Downward rounding")
    return dict(case=box["case"],block_radius_bits=bits,precision_bits=ctx.prec,iterates=N,
        native_source_radius_squared_max_approx=float(radius_squared.upper()),
        native_minimum_slab_gap_strictly_above=str(gap_lower),squared_shear_coefficients_exact=terms,
        one_inverse_Lipschitz_upper_approximation=float(L.upper()),
        return_inverse_Lipschitz_upper_approximation=float(invG.upper()),
        strict_rational_dimension_lower=str(lower),dimension_lower=float(lower))

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("out",type=Path)
    ap.add_argument("--bits",type=int,required=True);ap.add_argument("--source",type=Path,default=INPUT)
    a=ap.parse_args();r=bound(json.loads(a.source.read_text())["box"],a.bits)
    a.out.write_text(json.dumps(r,indent=2)+"\n");print(r["dimension_lower"],flush=True)

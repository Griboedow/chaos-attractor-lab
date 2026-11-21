"""Certify a dimension lower bound from backward binary itineraries on every slice."""
import argparse
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from verify_boxes import frames


def ball(q):
    return arb(q.numerator)/q.denominator


def bound(box,symbolic,trap,frame):
    bits=symbolic["parameter_radius_bits"]
    if trap["parameter_radius_bits"]!=bits or box["case"]!=trap["case"] or box["case"]!=symbolic["case"]:
        raise ValueError("The domain of the two computations must agree")
    if not symbolic.get("all_faces_and_full_slabs_verified"):
        raise ValueError("Complete two-slab crossing verification required")
    steps=trap["phase_regions"]
    if steps!=symbolic["cover"]["map_steps"] or steps!=box["case_order"]:
        raise ValueError("Different return maps")
    ctx.prec=128
    if frame["params"]!=box["params"] or frame["case"]!=box["case"]:
        raise ValueError("Slab and return maps differ")
    _,_,inverse=frames(frame)
    axis=symbolic["cover"]["axis"]
    gradient=sum(inverse[0][axis,j]**2 for j in range(3)).sqrt()
    faces=[arb(float(v)) for v in symbolic["cover"]["faces"]]
    separation=(faces[2]-faces[1])/gradient
    p=[Fraction.from_float(x) for x in box["params"]]
    radius=Fraction(0) if bits is None else Fraction(1,2**bits)
    rlo=p[0]-radius
    if not (0<rlo and p[0]+radius<1):
        raise ValueError("The proof requires affine contraction in the unit ball")
    factors=[]
    Lipschitz=1/ball(rlo)
    for i in range(3):
        k=abs(p[4+i])+radius;h=3*(abs(p[7+i])+radius)
        # Bound (k + h y) sqrt(1-y) on 0<=y<=1 for axis coordinate squared y.
        yy=[Fraction(0),Fraction(1)]
        if h:
            vertex=(2*h-k)/(3*h)
            if 0<vertex<1:
                yy.append(vertex)
        maximum=max((k+h*y)**2*(1-y) for y in yy)
        v=ball(maximum).sqrt()
        op=((4+v**2).sqrt()+v)/2
        factors.append(dict(max_shear_squared_exact=str(maximum),norm_upper_approximation=float(op.upper())))
        Lipschitz*=op
    returned=Lipschitz**steps
    lower=1+arb(2).log()/returned.log()
    rounded=Fraction(int(float(lower.lower())*1_000_000),1_000_000)
    if not lower>ball(rounded):
        raise ArithmeticError("Lower rounding failed")
    return dict(case=box["case"],parameter_radius_bits=bits,map_steps=steps,
        strict_native_slab_gap_bound_approximation=float(separation.lower()),
        longitudinal_gradient_norm_upper_approximation=float(gradient.upper()),
        first_slice_interval=symbolic["cover"]["faces"][:2],all_four_faces=symbolic["cover"]["faces"],
        one_step_inverse_Lipschitz_upper_approximation=float(Lipschitz.upper()),
        return_inverse_Lipschitz_upper_approximation=float(returned.upper()),
        rigorous_dimension_lower_rational=str(rounded),dimension_lower=float(rounded),
        reverse_twist_operator_bounds=factors,arb_precision_bits=ctx.prec)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for arg in ("box","symbolic","trap","frame"):
        ap.add_argument(arg,type=Path)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    result=bound(*(json.loads(getattr(a,k).read_text()) for k in ("box","symbolic","trap","frame")))
    a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(result["case"],result["dimension_lower"],flush=True)

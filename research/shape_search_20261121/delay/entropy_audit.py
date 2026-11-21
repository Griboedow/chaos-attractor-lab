"""Independently recheck rectangle coverage in separate x and y variables with Arb."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
from flint import arb,ctx

from search import save


def exact(q):
    p=Q(q)
    return arb(p.numerator)/p.denominator


def span(a,b):
    return exact(a).union(exact(b))


def check():
    ctx.prec=200
    r=span("4.48","4.52")
    a=span(".034",".035");b=span(".024",".025")
    c=1-a-b;s=span(".00045",".00055")
    low,high=Q(".0037"),Q("2.85")
    zlo,zhi=span(".0037",".209"),span(".38","2.85")
    margins=[1-r*c*zlo,r*c*zhi-1]
    cells=16
    for i in range(cells):
        x=span(low+(high-low)*i/cells,low+(high-low)*(i+1)/cells)
        for j in range(cells):
            y=span(low+(high-low)*j/cells,low+(high-low)*(j+1)/cells)
            def h(z):
                return s*x+(1-s)*z*(r*(1-a*x-b*y-c*z)).exp()
            margins.extend([exact(".38")-h(exact(".0037")),h(exact(".209"))-exact("2.85"),
                            h(exact(".38"))-exact("2.85"),exact(".0037")-h(exact("2.85"))])
    assert all(t>0 for t in margins)
    # Check the positive global minimum at each z endpoint using the worst
    # independent crowding coordinates. The survival contribution uses x_lo.
    floor=exact("1e-22")
    ceiling=exact(10)
    trapping=[]
    for z in [floor,ceiling]:
        g=(1-s)*z*(r*(1-a*ceiling-b*ceiling-c*z)).exp()
        trapping.append(exact(".00045")*floor+g-floor)
    trapping.append(ceiling-exact(".00055")*ceiling-(1-s)*(r-1).exp()/(r*c))
    assert all(t>0 for t in trapping)
    return {"precision_bits":ctx.prec,"independent_xy_boxes":cells**2,"strict_comparisons":len(margins)+len(trapping),
            "minimum_symbolic_margin_lower_approximation":float(min(x.lower() for x in margins)),
            "trapping_margins_lower_approximation":[float(x.lower()) for x in trapping]}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    result=check()
    save(args.out,result)
    print(result)

"""Bound the transverse multiplier of the exact closed ray by scalar Jacobi comparison."""
import argparse
import json
from pathlib import Path
from flint import arb, ctx


def certify(cells=16384):
    ctx.prec=128
    n=7
    k=arb(3)/8
    dt=arb.pi()*(arb(2)+arb(4)/n**2).sqrt()/cells
    low=(arb(1),arb(0),arb(0),arb(1))
    high=low

    def advance(vector, bound):
        u=bound.sqrt()*dt
        c=u.cosh()
        s=dt if bound==0 else dt*u.sinh()/u
        a,b,e,f=vector
        return (c*a+s*e,c*b+s*f,bound*s*a+c*e,bound*s*b+c*f)

    for i in range(cells):
        x=arb.pi()/2+arb.pi()*arb(2*i+1,1)/(2*cells)
        q=(n*k*x.cos()/(2*x.sin()**2+(n*n+2)*k*k))**2
        low=advance(low,max(arb(0),q.lower()))
        high=advance(high,q.upper())
    def full_ray(m):
        trace=(m[0]+m[3])/2
        return (trace+(trace**2-1).sqrt())**(2*n)
    lo,hi=full_ray(low),full_ray(high)
    assert lo>4450000 and hi<4470000
    return dict(cells=cells,bits=ctx.prec,half_wave_duration=str(dt*cells),
                lower_comparison_multiplier=str(lo),upper_comparison_multiplier=str(hi),
                certified_strict_bounds=[4450000,4470000])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    result=certify()
    args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result)

"""Enclose selected diagonal equilibria of cyclic harmonic flows."""
import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq


def binary_interval(value):
    numerator, denominator = float(value).as_integer_ratio()
    return mp.iv.mpf(numerator)/denominator


def interval_endpoints(value):
    # Preserve useful decimal bounds by rounding each side one extra ulp out.
    lo = np.nextafter(float(value.a), -np.inf)
    hi = np.nextafter(float(value.b), np.inf)
    return [float(lo), float(hi)]


def equilibria(case):
    p=np.asarray(case['parameters'])
    if p.shape != (3,9) or np.any(p!=p[0]):
        return []
    d,a,b,c,ph,ps,ch,m,n=p[0]
    def scalar(x):
        return -d*x+a*np.sin(x+ph)+b*np.sin(x+ps)+c*np.sin((m+n)*x+ch)
    bound=(abs(a)+abs(b)+abs(c))/d
    xs=np.linspace(-bound,bound,8193)
    ys=scalar(xs)
    centers=[brentq(scalar,xs[i],xs[i+1]) for i in np.flatnonzero(ys[:-1]*ys[1:]<0)]
    centers.extend(xs[np.abs(ys)<1e-13].tolist())
    centers=sorted(set(round(x,10) for x in centers))
    d,a,b,c,ph,ps,ch,m,n=map(binary_interval,p[0])
    def scalar_iv(x):
        return -d*x+a*mp.iv.sin(x+ph)+b*mp.iv.sin(x+ps)+c*mp.iv.sin((m+n)*x+ch)
    out=[]
    for center in centers:
        left, right = center-1e-9, center+1e-9
        lo, hi=binary_interval(left), binary_interval(right)
        f_lo, f_hi=scalar_iv(lo),scalar_iv(hi)
        if not (f_lo.a>0 and f_hi.b<0 or f_lo.b<0 and f_hi.a>0):
            continue
        region=mp.iv.mpf([left,right])
        gj=a*mp.iv.cos(region+ph)+m*c*mp.iv.cos((m+n)*region+ch)
        gk=b*mp.iv.cos(region+ps)+n*c*mp.iv.cos((m+n)*region+ch)
        real=-d+gj+gk
        pair=-d-(gj+gk)/2
        imag=mp.iv.sqrt(3)*(gj-gk)/2
        if not(real.a>0 or real.b<0):
            continue
        out.append(dict(diagonal_x_interval=[left,right],f_left=interval_endpoints(f_lo),
                        f_right=interval_endpoints(f_hi),real_axis_eigenvalue=interval_endpoints(real),
                        transverse_real_part=interval_endpoints(pair),transverse_imag_part=interval_endpoints(imag),
                        sink=bool(real.b<0 and pair.b<0)))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    mp.iv.dps=60
    selected=json.loads(args.manifest.read_text(encoding='utf-8-sig'))['cases']
    report={}
    for case in selected:
        if case['family']=='cartesian':
            roots=equilibria(case)
            if roots:
                report[str(case['case'])]=roots
    args.out.write_text(json.dumps(dict(precision_digits=mp.iv.dps,parameter_interpretation='Exact binary64 constants',
        scope='Each listed box has one diagonal zero by opposite signs and a uniformly nonzero derivative. Other zeros were not excluded.',
        cases=report),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()

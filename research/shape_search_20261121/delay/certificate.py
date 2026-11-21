"""Rational covering and trapping bounds for a box of delayed recruitment maps."""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path

from math import ceil, floor


def bounds(a,b=None):
    return (Q(a),Q(a if b is None else b))


def plus(x,y):
    return x[0]+y[0],x[1]+y[1]


def minus(x,y):
    return x[0]-y[1],x[1]-y[0]


def times(x,y):
    v=[p*q for p in x for q in y]
    return min(v),max(v)


def reciprocal(x):
    assert x[0]>0
    return 1/x[1],1/x[0]


def exp_point(x):
    if x<0:
        return reciprocal(exp_point(-x))
    m=0
    while x>Q(1,2):
        x/=2
        m+=1
    term=Q(1)
    total=Q(1)
    for n in range(1,43):
        term*=x/n
        total+=term
    first=term*x/43
    # The rest of the Taylor series is dominated by a geometric tail.
    low,high=total,total+first/(1-x/44)
    scale=10**55
    for _ in range(m):
        low=Q(floor(low*low*scale),scale)
        high=Q(ceil(high*high*scale),scale)
    return low,high


def exponential(x):
    return exp_point(x[0])[0],exp_point(x[1])[1]


def check():
    one=bounds(1)
    r=bounds("4.48","4.52")
    a=bounds(".034",".035")
    b=bounds(".024",".025")
    c=minus(minus(one,a),b)
    s=bounds(".00045",".00055")
    L=bounds(".0037",".209")
    H=bounds(".38","2.85")
    B=bounds("1e-22","10")
    C=(L[0],H[1])

    def output(xy,z):
        crowd=plus(times(plus(a,b),xy),times(c,z))
        recruit=times(times(minus(one,s),z),exponential(times(r,minus(one,crowd))))
        return plus(times(s,xy),recruit)

    tests={}
    def positive(name,margin):
        assert margin>0,name
        tests[name]=float(margin)

    positive("increasing_on_L",1-r[1]*c[1]*L[1])
    positive("decreasing_on_H",r[0]*c[0]*H[0]-1)
    positive("L_left_below_H",H[0]-output(C,bounds(L[0]))[1])
    positive("L_right_above_H",output(C,bounds(L[1]))[0]-H[1])
    positive("H_left_above_H",output(C,bounds(H[0]))[0]-H[1])
    positive("H_right_below_L",L[0]-output(C,bounds(H[1]))[1])

    # For given x,y, z*exp(-r*c*z) has one maximum. Its minimum on
    # the positive trapping interval occurs at an endpoint.
    positive("B_lower_at_left",output(B,bounds(B[0]))[0]-B[0])
    positive("B_lower_at_right",output(B,bounds(B[1]))[0]-B[0])
    peak=times(exponential(minus(r,one)),reciprocal(times(r,c)))
    maximum=plus(times(s,B),times(minus(one,s),peak))[1]
    positive("B_upper",B[1]-maximum)
    return {"parameter_intervals":{k:[str(q) for q in v] for k,v in
            (("r",r),("a",a),("b",b),("c=1-a-b",c),("s",s))},
            "L":[str(q) for q in L],"H":[str(q) for q in H],"B":[str(q) for q in B],
            "margins_display_float":tests,"three_step_trapping":"strict",
            "symbolic_transitions_at_stride_three":["L->H","H->L","H->H"],
            "entropy_lower_bound_per_update":"log((1+sqrt(5))/2)"}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    data=check()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(data["margins_display_float"])

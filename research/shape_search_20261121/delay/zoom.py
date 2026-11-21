"""Resolve the 39-fold period-doubling window on the diagonal of the trapping box."""
import argparse
import csv
import math
from pathlib import Path

import numpy as np
from numba import njit
from continuation import parameters


@njit(cache=True)
def sample(p,start,burn,retain,max_period):
    r,a,b,c,s=p
    x,y,z=start
    tail=np.empty(retain)
    for i in range(burn+retain):
        x,y,z=s*x+(1-s)*z*math.exp(r*(1-a*x-b*y-c*z)),x,y
        if i>=burn:
            tail[i-burn]=x
    cand=39
    while cand<=max_period:
        error=np.max(np.abs(tail[-(cand+200):]-tail[-(2*cand+200):-cand]))
        if error<3e-8:
            return np.array([x,y,z]),cand,error
        cand*=2
    return np.array([x,y,z]),0,error


def scan(out,start,end,n,burn,max_period):
    out.parent.mkdir(parents=True,exist_ok=True)
    grids=np.linspace(start,end,n,dtype=np.int64)
    assert len(set(grids))==n
    den=10**9
    with out.open("w",newline="",encoding="utf-8") as f:
        log=csv.writer(f)
        log.writerow(["direction","numerator","denominator","period_multiple_of_39","comparison_error","x","y","z"])
        for direction,indices in [("forward",grids),("reverse",grids[::-1])]:
            initial=np.array([.59,1.12,.47])
            counts={}
            for ix in indices:
                p=parameters(ix/den)
                initial,per,error=sample(p,initial,burn,2*max_period+200,max_period)
                counts[per]=counts.get(per,0)+1
                log.writerow([direction,ix,den,per,error,*initial])
            print(direction,counts,flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--first",type=int,default=978000000)
    ap.add_argument("--last",type=int,default=991000000)
    ap.add_argument("--count",type=int,default=4097)
    ap.add_argument("--burn",type=int,default=50000)
    ap.add_argument("--max-period",type=int,default=1248)
    args=ap.parse_args()
    scan(args.out,args.first,args.last,args.count,args.burn,args.max_period)

"""Test long transients and alternative initial conditions at certified periodic windows."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numba import njit

from confirm import native_spectrum
from search import save


@njit(cache=True)
def capture(p,z,orbit,sort,limit):
    r,a,b,c,s=p
    x,y,w=z
    n=len(orbit)
    matched=-1
    count=0
    for i in range(1,limit+1):
        x,y,w=s*x+(1-s)*w*math.exp(r*(1-a*x-b*y-c*w)),x,y
        if matched>=0:
            matched=(matched+1)%n
            if max(abs(orbit[matched,0]-x),abs(orbit[matched,1]-y),abs(orbit[matched,2]-w))<1e-7:
                count+=1
                if count>=3*n:
                    return i,np.array([x,y,w])
            else:
                matched=-1
        elif i%71==0:
            lo=np.searchsorted(sort[:,0],x-1e-8)
            while lo<n and sort[lo,0]<x+1e-8:
                k=int(sort[lo,1])
                if abs(orbit[k,1]-y)<1e-8 and abs(orbit[k,2]-w)<1e-8:
                    matched=k
                    count=1
                    break
                lo+=1
    return 0,np.array([x,y,w])


def probe(item,nstart,limit):
    t=item["t_numerator"]/item["t_denominator"]
    a=(34+t)/1000
    b=(24+t)/1000
    p=np.array([(448+4*t)/100,a,b,1-a-b,(45+10*t)/100000])
    r,a,b,c,s=p
    per=item["period"]
    x,y,z=map(float,item["center"])
    cycle=[]
    for _ in range(per):
        cycle.append((x,y,z))
        x,y,z=s*x+(1-s)*z*math.exp(r*(1-a*x-b*y-c*z)),x,y
    cycle=np.asarray(cycle)
    err=np.max(np.abs(cycle[0]-[x,y,z]))
    assert err<1e-8
    order=np.argsort(cycle[:,0])
    sorter=np.column_stack([cycle[order,0],order])
    draws=np.random.default_rng(2026112129).uniform(1e-8,10,(nstart,3))
    ts=[]
    active=[]
    for ic in draws:
        hit,last=capture(p,ic,cycle,sorter,limit)
        if hit:
            ts.append(hit)
        elif len(active)<12:
            # Continue the very same state with an independent tangent computation.
            rates,logs,_=native_spectrum(p,last,80000,8,100000,1)
            assert np.max(np.abs(rates.sum(1)-logs))<1e-6
            additional,_=capture(p,last,cycle,sorter,50000000)
            active.append({"initial":ic.tolist(),"at_limit":last.tolist(),"subsequent_spectra":rates.tolist(),
                           "extended_capture_time":limit+int(additional) if additional else None,
                           "extended_cutoff":limit+50000000})
    return {"period":per,"t":t,"draws":nstart,"limit":limit,"matches":len(ts),
            "double_precision_cycle_closure":float(err),
            "capture_time_quantiles_0_50_90_100":np.quantile(ts,[0,.5,.9,1]).tolist() if ts else None,
            "long_unmatched_examples":active}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificates",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--periods",type=int,nargs="+",default=[7,21,78,156])
    ap.add_argument("--count",type=int,default=256)
    ap.add_argument("--limit",type=int,default=1500000)
    args=ap.parse_args()
    rows=json.loads(args.certificates.read_text())
    results=[]
    for row in rows:
        if row["period"] in args.periods:
            rec=probe(row,args.count,args.limit)
            results.append(rec)
            save(args.out,results)
            print(rec["period"],rec["matches"],rec["draws"],flush=True)

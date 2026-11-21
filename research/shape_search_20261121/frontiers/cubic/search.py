"""Screen real Fourier drives of a cubic oscillator for emerging three-component links."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path

import numpy as np
from numba import njit


@njit(cache=True)
def mul(a,b):
    return a[0]+b[0],a[1]+b[1],a[2]+b[2]+a[0]*b[1]


@njit(cache=True)
def inv(a):
    return -a[0],-a[1],a[0]*a[1]-a[2]


@njit(cache=True)
def classify(values):
    """Track the three braid strands and the third longitude in the Heisenberg group."""
    labels=np.arange(3)
    counts=np.zeros((3,3),np.int64)
    x=[(1,0,0),(0,1,0),(0,0,0)]
    conjugators=[(0,0,0),(0,0,0),(0,0,0)]
    crossings=0
    for j,real_part in enumerate(values):
        if abs(real_part)<=2:
            continue
        slot=1 if real_part>2 else 0
        sign=(1 if j%2==0 else -1)*(1 if slot==1 else -1)
        i,k=labels[slot],labels[slot+1]
        counts[min(i,k),max(i,k)]+=sign
        p,q=conjugators[slot],conjugators[slot+1]
        if sign>0:
            m=mul(mul(inv(p),x[i]),p)
            conjugators[slot]=mul(q,inv(m))
            conjugators[slot+1]=p
        else:
            m=mul(mul(inv(q),x[k]),q)
            conjugators[slot]=q
            conjugators[slot+1]=mul(p,m)
        labels[slot],labels[slot+1]=k,i
        crossings+=1
    if not (labels[0]==0 and labels[1]==1 and labels[2]==2):
        return -1,0,crossings
    if counts[0,1] or counts[1,2] or counts[0,2]:
        return -2,0,crossings
    c=conjugators[2]
    assert c[0]==0 and c[1]==0
    return 0,c[2],crossings


@njit(cache=True)
def batch(values):
    out=np.zeros((len(values),3),np.int64)
    for i in range(len(values)):
        out[i]=classify(values[i])
    return out


def basis(t,h):
    angles=np.outer(t,np.arange(1,h+1))
    return np.column_stack((np.ones(len(t)),np.cos(angles),np.sin(angles)))


def search(n,seed):
    rng=np.random.default_rng(seed)
    counts=Counter()
    top={}
    draws=0
    for k in (3,4,5,6,7,8,9,10,12,14,16):
        for h in (2,3,4,5,6,8,10):
            theta=np.arange(2*k)*np.pi/k
            matrix=basis(theta,h)
            for offset in range(0,n,2048):
                count=min(2048,n-offset)
                c=rng.integers(-23,24,(count,2*h+1))
                c[:,0]//=3
                real=(c@matrix.T)/8
                outcome=batch(real)
                margins=np.min(abs(abs(real)-2),axis=1)
                counts["near_cut_threshold"]+=int(np.sum(margins<=1.e-9))
                for status,num in zip(*np.unique(outcome[margins>1.e-9,0],return_counts=True)):
                    counts[str(status)]+=int(num)
                winners=np.flatnonzero((outcome[:,0]==0)&(outcome[:,1]!=0)&(margins>.015))
                counts["nonzero_mu_screen_hits"]+=len(winners)
                for i in winners:
                    mu=int(outcome[i,1])
                    score=float(margins[i])/(1+.065*h+.015*k)
                    key=(abs(mu),h,k)
                    if key not in top or score>top[key]["score"]:
                        top[key]=dict(case=int(41000000+draws+i),h=h,k=k,coeffs=c[i].tolist(),divisor=8,
                                      imag_amplitude=3,mu=mu,crossings=int(outcome[i,2]),score=score,
                                      sampled_threshold_margin=float(margins[i]))
                draws+=count
    candidates=sorted(top.values(),key=lambda r:(abs(r["mu"]),r["h"],r["k"]))
    return dict(seed=seed,draws=draws,per_degree_pair=n,totals=dict(counts),best_by_degree_and_mu=candidates)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--per-pair",type=int,default=8192);ap.add_argument("--seed",type=int,default=220113)
    ap.add_argument("--out",type=Path,required=True);args=ap.parse_args()
    data=search(args.per_pair,args.seed)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(data["draws"],data["totals"], sorted({abs(r["mu"]) for r in data["best_by_degree_and_mu"]}),flush=True)

"""Adapt a return box to interior extrema located between regular-grid samples."""
import argparse
import copy
import json
from pathlib import Path
from box_scout import points_step
from cover_scout import charts
from validate import write_json
import numpy as np
from scipy.optimize import minimize


def extreme(box):
    o,M,iv=charts(box);p=np.array(box["params"]);m=box["case_order"]
    rng=np.random.default_rng(673352)
    g=np.linspace(-1,1,13)
    u=np.concatenate((rng.uniform(-1,1,(6000,3)),np.stack(np.meshgrid(g,g,g),axis=-1).reshape(-1,3)))
    x=o+u@M.T;avoid=np.inf
    for j in range(m):
        x=points_step(x,p)
        v=(x-o)@iv.T
        if j+1<m:
            avoid=min(avoid,float(np.min(np.max(abs(v),axis=1))))
    low=[];high=[]
    for ax in range(3):
        bounds=[]
        for sign in (1,-1):
            starts=[]
            for idx in np.argsort(sign*v[:,ax]):
                if all(np.max(abs(s-u[idx]))>.35 for s in starts):
                    starts.append(u[idx])
                    if len(starts)==10:
                        break
            def loss(z):
                point=(o+M@z).reshape(1,3)
                for _ in range(m):
                    point=points_step(point,p)
                return sign*float((iv@(point[0]-o))[ax])
            best=min(minimize(loss,z,method="L-BFGS-B",bounds=[(-1,1)]*3,
                              options=dict(ftol=1e-14,gtol=1e-10,maxiter=120)).fun for z in starts)
            bounds.append(sign*best)
        low.append(bounds[0]);high.append(bounds[1])
    return np.array(low),np.array(high),avoid


def expand(original):
    b=copy.deepcopy(original)
    source_lo=np.array(b["lo"][0]);source_hi=np.array(b["hi"][0])
    for k in range(36):
        lo=np.array(b["lo"][0]);hi=np.array(b["hi"][0]);mid=(lo+hi)/2;half=(hi-lo)/2
        low,high,avoid=extreme(b)
        margin=float(min(np.min(low+1),np.min(1-high)))
        if min(margin,avoid-1)>.008:
            b.update(optimizer_return_margin=margin,optimizer_avoidance_margin=avoid-1,optimizer_passes=k+1)
            return b
        if avoid<1.10 or np.max((hi-lo)/(source_hi-source_lo))>5:
            return None
        proposed_low=mid+half*(low-.035);proposed_high=mid+half*(high+.035)
        lo=np.where(low<-.985,np.minimum(lo,proposed_low),lo)
        hi=np.where(high>.985,np.maximum(hi,proposed_high),hi)
        b["lo"]=[lo.tolist()];b["hi"]=[hi.tolist()]
    return None


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("box",type=Path);ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args();b=json.loads(a.box.read_text());r=expand(b)
    if r:
        write_json(a.output,r)
        print(r["case"],r["optimizer_passes"],r["optimizer_return_margin"],flush=True)
    else:
        print(b["case"],"No optimized box",flush=True)

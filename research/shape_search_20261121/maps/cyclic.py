"""Search for positive expansion inside phase-separated cycles of spatial clouds."""
import argparse
import json
from pathlib import Path
from search import CACHE, START, design, sample
import numpy as np
from numba import njit


@njit(cache=True)
def scan_order(z,maximum):
    n=len(z)
    mean=z.sum(axis=0)/n
    centered=z-mean
    total=np.sum(centered*centered)
    if total < 1e-13:
        return 1,0.
    values=np.zeros(maximum+1)
    for m in range(2,maximum+1):
        sums=np.zeros((m,3))
        counts=np.zeros(m)
        for j in range(n):
            k=j%m
            counts[k]+=1
            for a in range(3):
                sums[k,a]+=centered[j,a]
        for k in range(m):
            for a in range(3):
                values[m]+=sums[k,a]*sums[k,a]/counts[k]/total
    best=values.max()
    for m in range(2,maximum+1):
        if values[m]>=best-.021:
            return m,values[m]
    return 1,0.


def phase_metrics(z,m):
    stop=len(z)-len(z)%m
    cycles=z[:stop].reshape(-1,m,3)
    centers=cycles.mean(axis=0)
    centered=z-z.mean(axis=0)
    residual=np.linalg.norm(cycles-centers,axis=2)
    radius=np.quantile(residual,.95,axis=0)
    delta=centers[:,None,:]-centers[None,:,:]
    distance=np.linalg.norm(delta,axis=2)
    distance.flat[::m+1]=np.inf
    ratios=distance/np.maximum(radius[:,None]+radius[None,:],1e-16)
    return dict(order=m,variance_explained=float(1-np.mean(residual**2)/np.mean(np.linalg.norm(centered,axis=1)**2)),
                minimum_center_95radius_ratio=float(ratios.min()),median_nearest_ratio=float(np.median(ratios.min(axis=1))),
                phase_radii_95=radius.tolist(),centers=centers.tolist())


def find_batch(cache,batch,count,ceiling):
    out=cache/f"cycle-{batch:03}"
    out.mkdir(parents=True,exist_ok=True)
    params=design(batch,count)
    hits=[]
    npositive=0
    for begin in range(0,count,128):
        p=params[begin:begin+128]
        ws,tails=sample(p,np.tile(START,(len(p),1)),3500,8192,8,2048)
        for j in np.where((ws.min(axis=1)>.008)&(ws.mean(axis=1)<.65))[0]:
            npositive+=1
            z=tails[j]
            m,v=scan_order(z,ceiling)
            if v<.60:
                continue
            info=phase_metrics(z,m)
            if info["median_nearest_ratio"]<.50:
                continue
            # Retain small clouds as well; the longer pass will characterize
            # high-order chains whose components may barely be visible at scale.
            hits.append(dict(case=1_000_000+batch*100_000+begin+int(j),params=p[j].tolist(),rates=ws[j].tolist(),**info))
    (out/"screen.json").write_text(json.dumps(dict(batch=batch,count=count,screen_positive=npositive,hits=hits))+"\n")
    print(batch,count,npositive,len(hits),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE)
    p.add_argument("--first",type=int,default=8)
    p.add_argument("--batches",type=int,default=16)
    p.add_argument("--count",type=int,default=4096)
    p.add_argument("--max-order",type=int,default=36)
    a=p.parse_args()
    for b in range(a.first,a.first+a.batches):
        find_batch(a.cache,b,a.count,a.max_order)

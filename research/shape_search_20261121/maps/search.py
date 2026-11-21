"""Geometry search in the repository's dissipative three-axis twist model."""
import argparse
import json
import os
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[key] = "1"
os.environ["NUMBA_NUM_THREADS"] = str(max(1,min(2,int(os.getenv("SHAPE_SEARCH_THREADS","2")))))
import numpy as np
from numba import njit, prange

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CACHE = Path("D:/Codex/chaos-attractor-cache/shape-search-20261121/maps")
RESULTS = REPO / "results/shape_search_20261121/maps"
START = np.array([.1,.2,.3])
FAMILY = "Dissipative three-axis twist map"


@njit(cache=True)
def advance(p, z, v):
    for j in range(3):
        z[j] *= p[0]
        v[j] *= p[0]
    z[0] += 1. - p[0]
    for i in range(3):
        j,k = (i+1)%3,(i+2)%3
        w = z[k]
        cubic = p[7+i]*w*w
        angle = p[1+i] + (p[4+i]+cubic)*w
        c,s = np.cos(angle),np.sin(angle)
        u,y = c*z[i]-s*z[j],s*z[i]+c*z[j]
        da = (p[4+i]+3.*cubic)*v[k]
        v[i],v[j] = c*v[i]-s*v[j]-y*da,s*v[i]+c*v[j]+u*da
        z[i],z[j] = u,y
    norm = np.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
    for i in range(3):
        v[i] /= norm
    return np.log(norm)


@njit(cache=True,parallel=True)
def sample(params, starts, burn, observed, windows, keep):
    rates = np.zeros((len(params),windows))
    tails = np.empty((len(params),keep,3))
    width = observed//windows
    for i in prange(len(params)):
        z,v = starts[i].copy(),np.array([1.,0.,0.])
        for n in range(burn+observed):
            rate = advance(params[i],z,v)
            if n >= burn:
                rates[i,(n-burn)//width] += rate/width
            if n >= burn+observed-keep:
                tails[i,n-burn-observed+keep] = z
    return rates,tails


def design(batch,count):
    rng = np.random.default_rng(np.random.SeedSequence([20261121,805,batch]))
    q = rng.uniform(-1,1,(count,10))
    p = q.copy()
    p[:,0] = .34 + .54*(q[:,0]+1)/2
    p[:,1:4] *= np.pi
    p[:,4:7] *= 15.
    p[:,7:10] *= 24.
    p[1::4,7:10] = 0.
    p[2::4,4:7] *= .25
    p[2::4,7:10] *= 1.7
    p[3::4,4:6] *= .12
    p[3::4,7:9] *= .18
    return p


def aligned(tail):
    x = tail-tail.mean(axis=0)
    vals,axes = np.linalg.eigh(x.T@x/len(tail))
    if vals[-1]<1e-13:
        raise ValueError("collapsed cloud")
    x = x @ axes[:,::-1]/np.sqrt(vals[-1])
    signs = np.sign((x**3).mean(axis=0))
    x *= np.where(signs==0,1,signs)
    return x,np.maximum(vals[::-1],0)/vals[-1]


def shape_features(tail):
    x,values = aligned(tail)
    radii = np.linalg.norm(x,axis=1)
    q = np.quantile(radii,[.05,.25,.50,.75,.95,.99])
    span = np.quantile(x,.99,axis=0)-np.quantile(x,.01,axis=0)
    f = [*values[1:],*q[:5]/max(q[-1],1e-10),
         *(np.abs(np.mean(x**3,axis=0))/np.maximum(values,1e-8)**1.5),
         span[1]/max(span[0],1e-10),span[2]/max(span[0],1e-10)]
    for grid in [8,20,48]:
        box = np.clip(((x+3.)*grid/6.).astype(int),0,grid-1)
        idx = box[:,0]+grid*(box[:,1]+grid*box[:,2])
        counts = np.bincount(idx)
        occupied = counts[counts>0]
        prob = occupied/len(x)
        f.extend([len(occupied)/len(x),np.exp(-prob@np.log(prob))/len(x)])
    return np.array(f)


def scout(batch,count,cache=CACHE):
    root = cache / f"batch-{batch:03}"
    root.mkdir(parents=True,exist_ok=True)
    params = design(batch,count)
    case = 1_000_000+batch*100_000+np.arange(count)
    features = np.full((count,18),np.nan)
    rates = np.empty((count,4))
    for begin in range(0,count,128):
        p = params[begin:begin+128]
        w,tails = sample(p,np.tile(START,(len(p),1)),2048,4096,4,640)
        rates[begin:begin+len(p)] = w
        for j in np.flatnonzero((w.min(axis=1)>.012)&(w.mean(axis=1)<.73)):
            try:
                features[begin+j] = shape_features(tails[j])
            except ValueError:
                pass
    np.savez_compressed(root/"scout.npz",case=case,params=params,rates=rates,features=features)
    print(json.dumps(dict(batch=batch,count=count,positive=int((rates.min(axis=1)>.012).sum()),
                          shaped=int(np.isfinite(features[:,0]).sum()))),flush=True)


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE)
    p.add_argument("--first",type=int,default=0)
    p.add_argument("--batches",type=int,default=1)
    p.add_argument("--count",type=int,default=4096)
    a=p.parse_args()
    for batch in range(a.first,a.first+a.batches):
        scout(batch,a.count,a.cache)

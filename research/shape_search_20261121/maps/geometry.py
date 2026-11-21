"""Rigid-pose comparison with 96 historical twist clouds and current finalists."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
from search import CACHE, RESULTS, START, sample, aligned
from validate import write_json
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import qmc

ARCHIVE_SHA256="0bc04d3aff12727ab628c59ec174e1abc84d36dac714245190fc7f98b30d064e"
HISTORICAL=[3469,1327,2612,808,67,2595,711,372,1075,2291,2650,1572,3699,62,815,3554,
            1487,2111,3555,1964,3800,272,2259,721,307,1203,1265,3203,2531,3254,1700,1994,
            4083,2882,2483,3884,2809,365,1262,3036,1002,4036,2082,716,450,1,3084,3068,
            176,3852,1426,3751,2934,451,1348,933,2875,3268,1692,2444,3579,820,2374,
            1891,869,1670,1773,44,3913,2639,3171,3400,578,3094,3386,2675,2419,140,
            1792,3537,3522,1164,54,868,2700,2104,653,1162,715,619,3720,1709,1980,3417,1717,2175]
AXES=[(0,1),(0,2),(1,2)]
SIGNS=list(itertools.product((1,-1),repeat=3))


def historical(cache,archive):
    if archive:
        if hashlib.sha256(archive.read_bytes()).hexdigest()!=ARCHIVE_SHA256:
            raise ValueError("Unexpected historical packet hash")
        with np.load(archive) as old:
            if not np.array_equal(old["ids"],HISTORICAL):
                raise ValueError("Historical case IDs do not match")
            return old["tail"],"original-archive"
    u=qmc.Sobol(10,scramble=True,seed=20261104).random_base2(12)
    u[:,0]=.35+.58*u[:,0];u[:,1:4]=2*np.pi*u[:,1:4]-np.pi
    u[:,4:7]=24*u[:,4:7]-12;u[:,7:]=16*u[:,7:]-8;u[:2048,7:]=0
    _,clouds=sample(u[HISTORICAL],np.tile(START,(len(HISTORICAL),1)),100_000,50_000,10,50_000)
    return clouds,"recomputed-historical-parameters"


def histogram(z):
    x,v=aligned(z.astype(np.float64))
    n=32
    pix=np.clip(((x+3)*n/6).astype(int),0,n-1)
    planes=[]
    for i,j in AXES:
        g=np.bincount(pix[:,i]*n+pix[:,j],minlength=n*n).reshape(n,n).astype(float)
        g=gaussian_filter(g,.65)
        planes.append(np.sqrt(g/g.sum()))
    return np.array(planes)


def distances(h):
    variants=np.array([[np.concatenate([item[j,::s[a],::s[b]].ravel()
                                        for j,(a,b) in enumerate(AXES)]) for item in h] for s in SIGNS])
    flat=variants[0]
    similarity=np.maximum.reduce([flat@other.T/3 for other in variants])
    return np.sqrt(np.maximum(1-similarity,0))


def compare(cache,results,archive):
    ref,mode=historical(cache,archive)
    checks=[json.loads(p.read_text()) for p in sorted((results/"checks").glob("*.json"))]
    checks=[p for p in checks if p.get("cloud")]
    first=[histogram(z[:25_000]) for z in ref]
    second=[histogram(z[-25_000:]) for z in ref]
    for row in checks:
        with np.load(row["cloud"]) as cloud:
            z=cloud["points"]
            first.append(histogram(z[:25_000]))
            second.append(histogram(z[-25_000:]))
    n=len(HISTORICAL)
    all_dist=distances(np.array(first+second))
    d=all_dist[:len(first),:len(first)].copy()
    within=np.diag(all_dist[:len(first),len(first):])
    np.fill_diagonal(d,np.inf)
    historical_nn=d[:n,:n].min(axis=1)
    ids=HISTORICAL + [c["case"] for c in checks]
    comparisons=[]
    for i,c in enumerate(checks,n):
        j=int(np.argmin(d[i,:n]))
        k=n+int(np.argmin(d[i,n:])) if len(checks)>1 else None
        comparisons.append(dict(case=c["case"],nearest_archived_case=ids[j],archive_distance=float(d[i,j]),
            archived_leave_one_out_percentile=float(np.mean(historical_nn<=d[i,j])),
            own_disjoint_tail_distance=float(within[i]),
            nearest_current_case=(ids[k] if k is not None else None),
            current_distance=(float(d[i,k]) if k is not None else None)))
    packet=dict(archive_mode=mode,archive_sha256=ARCHIVE_SHA256,historical_ids=HISTORICAL,
                metric="Minimum over PCA-axis sign flips of RMS Hellinger distance across three 32-bin orthogonal probability projections. PCA centers each cloud, rescales all axes by first-axis RMS and preserves aspect ratios. A Gaussian 0.65-bin blur reduces sampling noise. Near-degenerate eigenvalues and asymmetric density can affect comparison.",
                historical_nearest_distance_quantiles=np.quantile(historical_nn,[0,.25,.5,.75,1]).tolist(),
                historical_disjoint_split_quantiles=np.quantile(within[:n],[0,.25,.5,.75,1]).tolist(),
                comparisons=comparisons)
    write_json(results/"geometry.json",packet)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE)
    p.add_argument("--results",type=Path,default=RESULTS)
    p.add_argument("--archive",type=Path)
    a=p.parse_args()
    compare(a.cache,a.results,a.archive)

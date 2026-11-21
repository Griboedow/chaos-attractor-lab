"""Find candidate parallelepiped trapping cycles by sampled closure iteration.

This construction is exploratory. Only a subsequent outward-rounded verifier
can turn its sampled forward-inclusion tests into a certificate.
"""
import argparse
import itertools
import json
from pathlib import Path
from search import CACHE, START, sample
from validate import write_json
import numpy as np
from numba import njit


@njit(cache=True)
def points_step(points,p):
    result=points.copy()
    for n in range(len(result)):
        x=result[n]
        for i in range(3):
            x[i]*=p[0]
        x[0]+=1.-p[0]
        for i in range(3):
            j,k=(i+1)%3,(i+2)%3
            w=x[k]
            a=p[1+i]+p[4+i]*w+p[7+i]*w**3
            c,s=np.cos(a),np.sin(a)
            x[i],x[j]=c*x[i]-s*x[j],s*x[i]+c*x[j]
    return result


def frames(z,m):
    stop=len(z)-len(z)%m
    groups=z[:stop].reshape(-1,m,3).transpose(1,0,2)
    center=groups.mean(axis=1)
    means=groups-center[:,None,:]
    covariance=means.transpose(0,2,1)@means
    _,basis=np.linalg.eigh(covariance)
    projected=means@basis
    low,high=projected.min(axis=1),projected.max(axis=1)
    center+=np.einsum("mij,mj->mi",basis,(low+high)/2)
    return center,basis,np.maximum((high-low)/2,1e-13)


def try_boxes(p,z,m,dimension_floor,grid_points=5,max_steps=240):
    c,axes,h=frames(z,m)
    scale=h.max(axis=1,keepdims=True)
    h=np.maximum(h,dimension_floor*scale)
    lo,hi=-1.08*h.copy(),1.08*h.copy()
    source_lo,source_hi=lo.copy(),hi.copy()
    g=np.array(list(itertools.product(np.linspace(-1,1,grid_points),repeat=3)))
    after=np.roll(np.arange(m),-1)
    for k in range(max_steps):
        mid=(hi+lo)/2
        half=(hi-lo)/2
        x=(mid[:,None,:]+g[None,:,:]*half[:,None,:])@axes.transpose(0,2,1)+c[:,None,:]
        y=points_step(x.reshape(-1,3),p).reshape(x.shape)
        mapped=(y-c[after,None,:])@axes[after]
        low=np.roll(mapped.min(axis=1),1,axis=0)
        high=np.roll(mapped.max(axis=1),1,axis=0)
        clearance=float(np.min(np.minimum(low-lo,hi-high)/half))
        if clearance>.003:
            return dict(case_order=m,center=c.tolist(),basis=axes.tolist(),lo=lo.tolist(),hi=hi.tolist(),
                        sampled_clearance=clearance,passes=k+1,dimension_floor=dimension_floor)
        active=(low<lo+.008*half)|(high>hi-.008*half)
        lo=np.where(active,np.minimum(lo,low-.024*half),lo)
        hi=np.where(active,np.maximum(hi,high+.024*half),hi)
        if np.max((hi-lo)/(source_hi-source_lo))>200 or np.max(np.abs(hi-lo))>1.8:
            break
    return None


def main(cache,preview,limit):
    rows=json.loads(preview.read_text())
    rows.sort(key=lambda r:-r["minimum_box_gap_fraction"])
    out=cache/"box-pilot"
    out.mkdir(parents=True,exist_ok=True)
    found=0
    for row in rows[:limit]:
        m=row["order"];p=np.array(row["params"])
        _,z=sample(p[None,:],START[None,:],80_000,30_000,3,30_000)
        for floor in [.06,.25,.6]:
            box=try_boxes(p,z[0],m,floor)
            if box:
                box.update(case=row["case"],params=p.tolist())
                write_json(out/f"{row['case']}.json",box)
                print(row["case"],m,box["sampled_clearance"],flush=True)
                found+=1
                break
    print(dict(tested=min(limit,len(rows)),sample_closures=found),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE)
    p.add_argument("--preview",type=Path,required=True)
    p.add_argument("--limit",type=int,default=300)
    a=p.parse_args()
    main(a.cache,a.preview,a.limit)

"""Find a trapping candidate for one return without boxing intermediate phases."""
import argparse
import itertools
import json
from pathlib import Path
from box_scout import frames, points_step
from search import CACHE, START, sample
from validate import write_json
import numpy as np


def attempt(p,m,c,basis,h,floor,grid):
    h=np.maximum(h,floor*h.max())
    lo,hi=-1.10*h.copy(),1.10*h.copy()
    initial_width=hi-lo
    for k in range(90):
        mid,half=(hi+lo)/2,(hi-lo)/2
        x=(mid+grid*half)@basis.T+c
        clearance=np.inf
        for j in range(1,m+1):
            x=points_step(x,p)
            projected=((x-c)@basis-mid)/half
            if j<m:
                clearance=min(clearance,float(np.min(np.max(np.abs(projected),axis=1)-1)))
                if clearance<.05:
                    return None
        low,high=projected.min(axis=0),projected.max(axis=0)
        slack=float(np.min(np.minimum(low+1,1-high)))
        if slack>.012:
            return dict(mode="return",case_order=m,center=[c.tolist()],basis=[basis.tolist()],
                        lo=[lo.tolist()],hi=[hi.tolist()],sampled_return_margin=slack,
                        sampled_avoidance_margin=clearance,passes=k+1,dimension_floor=floor)
        amin=mid+half*low;amax=mid+half*high
        active=(low<-1+.025)|(high>1-.025)
        lo=np.where(active,np.minimum(lo,amin-.07*half),lo)
        hi=np.where(active,np.maximum(hi,amax+.07*half),hi)
        if np.max((hi-lo)/initial_width)>100 or np.max(hi-lo)>1.2:
            break
    return None


def main(cache,preview,limit,max_order):
    rows=[r for r in json.loads(preview.read_text()) if r["order"]<=max_order]
    rows.sort(key=lambda r:(r["order"],-r["minimum_box_gap_fraction"]))
    grid=np.array(list(itertools.product(np.linspace(-1,1,5),repeat=3)))
    out=cache/"return-pilot";out.mkdir(parents=True,exist_ok=True)
    hits=0
    for row in rows[:limit]:
        m=row["order"];p=np.array(row["params"])
        _,z=sample(p[None,:],START[None,:],80_000,30_000,3,30_000)
        centers,bases,widths=frames(z[0],m)
        found=False
        for j in np.argsort(widths.max(axis=1)):
            for floor in (.05,.2,.5,.01,1.):
                box=attempt(p,m,centers[j],bases[j],widths[j],floor,grid)
                if box:
                    box.update(case=row["case"],params=p.tolist(),source_phase=int(j))
                    write_json(out/f"{row['case']}.json",box)
                    print(row["case"],m,round(box["sampled_return_margin"],4),round(box["sampled_avoidance_margin"],4),flush=True)
                    hits+=1;found=True;break
            if found:
                break
    print(dict(tested=min(len(rows),limit),proposals=hits),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE)
    p.add_argument("--preview",type=Path,required=True)
    p.add_argument("--max-order",type=int,default=14)
    p.add_argument("--limit",type=int,default=300)
    a=p.parse_args()
    main(a.cache,a.preview,a.limit,a.max_order)

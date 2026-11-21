"""Search return-box faces for two disjoint full-crossing intervals."""
import argparse
import json
from pathlib import Path
from box_scout import points_step
from search import CACHE
from validate import write_json
import numpy as np


def charts(packet):
    c=np.array(packet["center"][0]);b=np.array(packet["basis"][0])
    lo,hi=np.array(packet["lo"][0]),np.array(packet["hi"][0])
    o=c+b@((lo+hi)/2)
    M=b*((hi-lo)/2)
    return o,M,np.linalg.inv(M)


def best_four(u,low,high):
    best=.006;choice=None
    n=len(u)
    for i in range(n-6):
        for h in range(i+6,n):
            if u[h]-u[i]<.04:
                continue
            below=u[i]-high
            above=low-u[h]
            li=int(above[i]>below[i]);rh=int(above[h]>below[h])
            ends=min(max(above[i],below[i]),max(above[h],below[h]))
            if ends<=best:
                continue
            # Each branch maps its endpoints to opposite sides of both slabs.
            left=below if li else above
            right=below if rh else above
            ids=np.arange(i+1,h-2)
            nextbest=np.maximum.accumulate(right[i+3:h][::-1])[::-1]
            scores=np.minimum(np.minimum(left[ids],nextbest),ends)
            pos=int(np.argmax(scores))
            if scores[pos]>best:
                j=int(ids[pos]);k=j+2+int(np.argmax(right[j+2:h]))
                best=float(scores[pos])
                choice=dict(faces=list(map(float,u[[i,j,k,h]])),side=[li,1-li,1-rh,rh],
                            sampled_minimum_margin=best)
    return choice


def find(packet,n=193,transverse=5,returns=8):
    o,basis,inv=charts(packet)
    p=np.array(packet["params"]);m=packet["case_order"]
    out=[]
    for axis in range(3):
        other=[i for i in range(3) if i!=axis]
        u=np.linspace(-1.,1.,n)
        g=np.linspace(-1.,1.,transverse)
        mesh=np.stack(np.meshgrid(u,g,g,indexing="ij"),axis=-1).reshape(-1,3)
        s=mesh[:,np.argsort([axis,*other])]
        x=o+s@basis.T
        for k in range(1,returns+1):
            for _ in range(m):
                x=points_step(x,p)
            h=((x-o)@inv.T)[:,axis].reshape(n,-1)
            hit=best_four(u,h.min(axis=1),h.max(axis=1))
            if hit:
                hit.update(case=packet["case"],axis=axis,returns=k,map_steps=m*k)
                out.append(hit)
    return sorted(out,key=lambda x:(x["returns"],-x["sampled_minimum_margin"]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("boxes",type=Path,nargs="+")
    p.add_argument("--cache",type=Path,default=CACHE)
    p.add_argument("--returns",type=int,default=8)
    a=p.parse_args()
    for path in a.boxes:
        packet=json.loads(path.read_text())
        hits=find(packet,returns=a.returns)
        write_json(a.cache/"cover-pilot"/f"{packet['case']}.json",hits)
        print(packet["case"],hits[:2],flush=True)

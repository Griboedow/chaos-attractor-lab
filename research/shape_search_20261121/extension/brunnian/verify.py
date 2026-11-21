"""Independently verify the pure-braid longitudes, sampled separation and data receipts."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from sample import centers

# Integer Heisenberg matrices have their three strict upper entries stored as (x,y,z).
# Here x is entry 12, y is entry 23 and z is entry 13.
def mul(a,b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2]+a[0]*b[1])


def inv(a):
    return (-a[0],-a[1],a[0]*a[1]-a[2])


def full_twists(word):
    ident=(0,0,0)
    starts=[(1,0,0),(0,1,0),ident]
    meridians=starts.copy()
    left=[ident]*3
    for block in word:
        i=abs(block)-1
        product=mul(meridians[i],meridians[i+1])
        transform=product if block>0 else inv(product)
        for j in (i,i+1):
            left[j]=mul(transform,left[j])
            meridians[j]=mul(mul(left[j],starts[j]),inv(left[j]))
    result=inv(left[2])
    assert result[:2]==(0,0)
    return result[2]


def check(designs,records,cache):
    rows=[]
    for design,record in zip(designs,records,strict=True):
        word=design["word"]
        mu=full_twists(word)
        assert mu==design["mu_123"]
        assert [word.count(k)-word.count(-k) for k in (1,2)]==[0,0]
        theta=2*np.pi*np.arange(16384)/16384
        xy=[np.stack(centers(theta,word,i),axis=1) for i in range(3)]
        gap=min(np.linalg.norm(xy[i]-xy[j],axis=1).min() for i,j in ((0,1),(1,2),(0,2)))
        assert gap>=1-2e-14
        file=cache/record["cloud_relative"]
        assert hashlib.sha256(file.read_bytes()).hexdigest()==record["sha256"]
        with np.load(file) as c:
            z,labels=c["points"],c["component"]
        assert np.isfinite(z).all() and list(z.shape)==record["shape"] and record["case"]==design["case"]
        for strand in range(3):
            subset=z[labels==strand]
            assert len(subset)==record["steps_per_component"]
            angles=np.arctan2(subset[:,1],subset[:,0])
            u,v=centers(angles,word,strand)
            offset=np.column_stack((np.hypot(subset[:,0],subset[:,1])-design["radius"]-u,subset[:,2]-v))
            assert np.linalg.norm(offset,axis=1).max()<design["width"]
        rows.append(dict(case=design["case"],heisenberg_mu=mu,minimum_sampled_center_gap=float(gap),
                         points=len(z),file_sha256=record["sha256"]))
    return rows


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--packet",required=True,type=Path);ap.add_argument("--cache",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    records=json.loads((args.packet/"clouds.json").read_text())
    cases=json.loads((args.packet/"designs.json").read_text())["cases"]
    checks=check(cases,records,args.cache)
    args.out.write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
    print(f"{len(checks)} independent integral longitude and cloud checks")

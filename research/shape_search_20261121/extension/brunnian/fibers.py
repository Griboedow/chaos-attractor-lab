"""Enumerate inverse angular histories and check the predicted separation of fiber disks."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def pull_back(angles,d,ripple,shift):
    center=angles-shift
    low=center-ripple
    high=center+ripple
    for _ in range(54):
        mid=.5*(low+high)
        lower=mid+ripple*np.sin(mid)<center
        low=np.where(lower,mid,low)
        high=np.where(lower,high,mid)
    lift=.5*(low+high)
    offsets=2*np.pi*np.arange(d)
    past=(lift[:,None]+offsets[None,:])/d
    return (past%(2*np.pi)).ravel()


def run(row,depth):
    rho=row["fiber_contraction"]
    b=row["fiber_offset"]
    report=[]
    for d in row["degrees"]:
        gap=2*b*math.sin(math.pi/d)-2*rho
        samples=[]
        for target in (.277,1.938,5.218):
            sums=np.zeros(1,dtype=complex)
            theta=np.asarray([target])
            levels=[]
            for k in range(depth):
                theta=pull_back(theta,d,row["base_ripple"],row["base_shift"])
                sums=np.repeat(sums,d)+(rho**k)*b*np.exp(1j*theta)
                coords=np.column_stack((sums.real,sums.imag))
                minimum=float(cKDTree(coords).query(coords,k=2)[0][:,1].min())
                assert minimum>=gap*rho**k and np.max(np.abs(sums))+rho**(k+1)<1
                levels.append(dict(depth=k+1,disks=len(theta),min_center_distance=minimum,
                                   analytic_min_center_distance=gap*rho**k,radius=rho**(k+1)))
            samples.append(dict(target=target,levels=levels))
        report.append(dict(degree=d,fiber_dimension=math.log(d)/math.log(1/rho),
                           attractor_dimension=1+math.log(d)/math.log(1/rho),
                           uniform_first_level_gap=gap,sampled_sections=samples))
    return report


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("designs",type=Path);ap.add_argument("--depth",type=int,default=7)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    row=json.loads(args.designs.read_text())["cases"][0]
    result=run(row,args.depth)
    args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print([(r["degree"],r["attractor_dimension"]) for r in result])

"""Measure finite capture times near the certified period-39 sink."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from numba import njit

from search import save


@njit(cache=True)
def fate(initial,cycle,limit):
    r,a,b,c,s=4.52,.035,.025,.94,.00055
    x,y,z=initial
    candidate=-1
    first=-1
    matched=0
    for i in range(1,limit+1):
        x,y,z=s*x+(1-s)*z*math.exp(r*(1-a*x-b*y-c*z)),x,y
        if candidate>=0:
            candidate=(candidate+1)%39
            if max(abs(x-cycle[candidate,0]),abs(y-cycle[candidate,1]),abs(z-cycle[candidate,2]))<1e-7:
                matched+=1
                if matched>=117:
                    return first,i,x,y,z
            else:
                candidate=-1
                matched=0
        elif i%43==0:
            for k in range(39):
                if max(abs(x-cycle[k,0]),abs(y-cycle[k,1]),abs(z-cycle[k,2]))<1e-8:
                    candidate,first,matched=k,i,1
                    break
    return -1,limit,x,y,z


def run(cert,out,count,limit):
    data=json.loads(cert.read_text())
    r,a,b,c,s=map(float,data["parameters"])
    x,y,z=map(float,data["center"])
    cycle=[]
    for i in range(39):
        cycle.append((x,y,z))
        x,y,z=s*x+(1-s)*z*math.exp(r*(1-a*x-b*y-c*z)),x,y
    cycle=np.asarray(cycle)
    assert max(abs(cycle[0]-[x,y,z]))<1e-10
    rng=np.random.default_rng(2026112113)
    summaries={}
    out.mkdir(parents=True,exist_ok=True)
    with (out/"basins.csv").open("w",newline="",encoding="utf-8") as f:
        col=csv.writer(f)
        col.writerow(["cohort","x","y","z","first_matching_check","stop","last_x","last_y","last_z"])
        for name in ("uniform_cube","log_uniform","symbol_intervals"):
            if name=="uniform_cube":
                starts=rng.uniform(1e-12,10,(count,3))
            elif name=="log_uniform":
                starts=np.exp(rng.uniform(math.log(1e-12),math.log(10),(count,3)))
            else:
                modes=rng.integers(0,2,(count,3))
                starts=np.where(modes,rng.uniform(.38,2.85,(count,3)),rng.uniform(.0037,.209,(count,3)))
            times=[]
            for initial in starts:
                hit,stop,x,y,z=fate(initial,cycle,limit)
                col.writerow([name,*initial,hit,stop,x,y,z])
                if hit>=0:
                    times.append(hit)
            summaries[name]={"draws":count,"observed_captures":len(times),"limit":limit,
                             "first_matching_check_quantiles":np.quantile(times,[0,.5,.9,.99,1]).tolist() if times else None}
    save(out/"summary.json",{"seed":2026112113,"initials_per_group":count,
                             "time_quantiles":[0,.5,.9,.99,1],"cohorts":summaries,
                             "test":"Near a reference phase at a 43-step checkpoint within 1e-8, then within 1e-7 for 117 consecutive steps in phase order."})
    print(summaries,flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificate",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--count",type=int,default=4096)
    ap.add_argument("--limit",type=int,default=150000)
    args=ap.parse_args()
    run(args.certificate,args.out,args.count,args.limit)

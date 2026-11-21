"""Look for stable periodic basins alongside the selected numerical or certified sets."""
import argparse
import json
from pathlib import Path
from search import CACHE,RESULTS,sample
from validate import write_json
import numpy as np


def screen(case,count):
    row=json.loads((RESULTS/"checks"/f"{case}.json").read_text())
    p=np.array(row["params"])
    rng=np.random.default_rng(np.random.SeedSequence([20261121,3317,case]))
    starts=rng.uniform(-1,1,(count,3))
    starts*=10**rng.uniform(-1,1.15,(count,1))
    starts[:16]=[r["start"] for r in row["basin"]]
    hits=[];negative=0
    for begin in range(0,count,128):
        x=starts[begin:begin+128]
        w,zz=sample(np.tile(p,(len(x),1)),x,8000,8192,4,1024)
        for i in np.where(w.max(axis=1)<-.0005)[0]:
            z=zz[i]
            periods=np.flatnonzero(np.linalg.norm(z[-1]-z[-513:-1][::-1],axis=1)<1e-7)+1
            if not len(periods):
                continue
            m=int(periods[0])
            if np.max(abs(z[-2*m:-m]-z[-m:]))>1e-7:
                continue
            negative+=1
            orbit=z[-m:]
            if any(r["period"]==m and np.min(np.linalg.norm(np.array(r["orbit"])-orbit[0],axis=1))<1e-5 for r in hits):
                continue
            hits.append(dict(period=m,orbit=orbit.tolist(),leading_windows=w[i].tolist(),start=x[i].tolist()))
    out=CACHE/"basins"/f"{case}.json"
    write_json(out,dict(case=case,initial_states=count,negative_periodic_states=negative,hits=hits))
    print(case,negative,[(h["period"],h["leading_windows"][0]) for h in hits],flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases",type=int,nargs="+",required=True)
    ap.add_argument("--count",type=int,default=4096)
    a=ap.parse_args()
    for case in a.cases:
        screen(case,a.count)

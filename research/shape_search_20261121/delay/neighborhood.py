"""Search the full spectrum near four geometrically selected recruitment-delay regimes."""
import argparse
import csv
import math
import json
from pathlib import Path

import numpy as np

from confirm import native_spectrum
from search import save


def run(screen,out,count):
    previous=json.loads(screen.read_text())["diversity_sample"]
    seeds=[5346,28316,30498,13753]
    anchors={v["case"]:np.array(v["parameters"]) for v in previous if v["case"] in seeds}
    assert set(anchors)==set(seeds)
    rng=np.random.default_rng(2026112110)
    initial=np.array([.59,1.12,.47])
    retained=[]
    out.mkdir(parents=True,exist_ok=True)
    with (out/"local-screen.csv").open("w",newline="",encoding="utf-8") as f:
        log=csv.writer(f)
        log.writerow(["case","anchor","r","a","b","c","s","min_l2","min_l3","det_error","cells"])
        for index in range(count):
            anchor=seeds[index%4]
            base=anchors[anchor]
            width=rng.uniform(.10,1.05)
            weights=base[1:4]*np.exp(rng.normal(0,width,3))
            weights/=weights.sum()
            r=np.clip(base[0]+rng.normal(0,.57),2.55,5.3)
            s=np.clip(base[4]*math.exp(rng.normal(0,.72)),.025,.52)
            par=np.r_[r,weights,s]
            ident=1000000+index
            try:
                rates,logs,z=native_spectrum(par,initial,2500,3,800,1200)
            except ValueError:
                continue
            det_error=float(np.max(np.abs(rates.sum(1)-logs)))
            if det_error>1e-6 or not np.isfinite(rates).all():
                continue
            l2,l3=rates.min(0)[1:]
            if l2<.008:
                continue
            scale=z.max()-z.min()
            if scale<1e-6:
                continue
            cells=np.floor((z-z.min())/scale*18).astype(np.int32)
            _,counts=np.unique(cells[:,0]*361+cells[:,1]*19+cells[:,2],return_counts=True)
            log.writerow([ident,anchor,*par,l2,l3,det_error,len(counts)])
            if len(counts)>80 and counts.max()<.14*len(z):
                retained.append({"case":ident,"parameters":par.tolist(),"anchor":anchor,
                                 "short_min_l2":float(l2),"short_min_l3":float(l3),"cells":len(counts)})
    # Preserve both high growth and moderate growth with strong coupling.
    keep=[]
    for score in (lambda r:r["short_min_l2"],lambda r:r["short_min_l3"],
                  lambda r:r["parameters"][4],lambda r:-r["cells"],
                  lambda r:min(r["parameters"][1:4])):
        for anchor in seeds:
            take=sorted((v for v in retained if v["anchor"]==anchor),key=score,reverse=True)[:6]
            keep.extend(v for v in take if v not in keep)
    save(out/"screen.json",{"seed":2026112110,"draws":count,"recruitment_delay":True,
                           "passed_occupancy":len(retained),"diversity_sample":keep})
    print(count,len(retained),len(keep))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("screen",type=Path)
    ap.add_argument("--out",required=True,type=Path)
    ap.add_argument("--count",default=24000,type=int)
    args=ap.parse_args()
    run(args.screen,args.out,args.count)

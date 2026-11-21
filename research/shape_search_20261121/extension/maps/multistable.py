"""Screen reproducibly perturbed neighborhoods of three multistable reference maps."""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/"research/shape_search_20261121/maps"))
import numpy as np
from search import START,sample

OLD=ROOT/"results/shape_search_20261121/maps"
CASES=(20100964,2603621,2303546)

def screen(cache,count):
    cache.mkdir(parents=True,exist_ok=True)
    basins={b["case"]:b for b in json.loads((OLD/"additional-basins.json").read_text())["cases"]}
    hits=[]
    for case in CASES:
        base=np.array(json.loads((OLD/"checks"/f"{case}.json").read_text())["params"])
        rng=np.random.default_rng(np.random.SeedSequence([20261121,2110,case]))
        param=np.tile(base,(count,1))
        direction=rng.normal(size=(count,10))
        direction*=np.array([.15,1,1,1,5,5,5,8,8,8])
        param+=direction*10**rng.uniform(-4,-1.35,(count,1))
        np.savez_compressed(cache/f"design-{case}.npz",params=param)
        for head in range(0,count,64):
            x=param[head:head+64];w,z=sample(x,np.tile(START,(len(x),1)),6000,4096,4,1024)
            for k in np.flatnonzero(w.mean(axis=1)<-.001):
                t=z[k];ret=np.flatnonzero(np.linalg.norm(t[-1]-t[-513:-1][::-1],axis=1)<1e-8)+1
                if len(ret)==0:continue
                m=int(ret[0])
                if m<8 or np.max(np.abs(t[-2*m:-m]-t[-m:]))>1e-8:continue
                hits.append(dict(id=f"{case}-{head+k:04d}",parent=case,index=head+int(k),period=m,
                    rate=float(w[k].mean()),params=x[k].tolist(),start=START.tolist(),
                    alternate_starts=basins[case]["example_starts"]))
            (cache/"scout.json").write_text(json.dumps(hits))
        print(case,count,len(hits),sorted({h["period"] for h in hits if h["parent"]==case}),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cache",type=Path);ap.add_argument("--count",type=int,default=384)
    a=ap.parse_args();screen(a.cache,a.count)

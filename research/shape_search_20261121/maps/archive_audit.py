"""Compare against every noncollapsed cloud in the original 4,096-case scout."""
import argparse
import hashlib
import json
from pathlib import Path
from search import CACHE, RESULTS, START, sample
from geometry import histogram, SIGNS, AXES
from cyclic import scan_order,phase_metrics
from validate import write_json
import numpy as np
from scipy.stats import qmc


SHA="595163d76571184092e1a37b6710dbc689f69a7b9e2bf427494db7f1e11b20d5"


def source(archive):
    if archive:
        if hashlib.sha256(archive.read_bytes()).hexdigest()!=SHA:
            raise ValueError("Historical scout digest mismatch")
        with np.load(archive) as data:
            return data["tail"],data["windows"],"original-archive"
    p=qmc.Sobol(10,scramble=True,seed=20261104).random_base2(12)
    p[:,0]=.35+.58*p[:,0];p[:,1:4]=2*np.pi*p[:,1:4]-np.pi
    p[:,4:7]=24*p[:,4:7]-12;p[:,7:]=16*p[:,7:]-8;p[:2048,7:]=0
    rates,tails=sample(p,np.tile(START,(len(p),1)),5000,12000,4,2048)
    return tails,rates,"recomputed-historical-scout"


def run(results,archive):
    tails,win,mode=source(archive)
    hs=[];ids=[];phases=[]
    for i,z in enumerate(tails):
        try:
            h=histogram(z)
        except ValueError:
            continue
        hs.append(h.astype("float32"));ids.append(i)
        if win[i].min()>.008:
            m,v=scan_order(z,36)
            if v>=.60 and phase_metrics(z,m)["median_nearest_ratio"]>=.5:
                phases.append(dict(case=i,order=m,phase_variance=float(v)))
    hist=np.array(hs).reshape(len(hs),-1)
    rows=[]
    for f in sorted((results/"checks").glob("*.json")):
        b=json.loads(f.read_text())
        if not b.get("cloud"):
            continue
        z=np.load(b["cloud"])["points"]
        h=histogram(z[:tails.shape[1]])
        signs=np.array([np.concatenate([h[j,::s[a],::s[c]].ravel() for j,(a,c) in enumerate(AXES)]) for s in SIGNS],dtype="float32")
        sim=np.max(hist@signs.T,axis=1)/3
        j=int(np.argmax(sim))
        rows.append(dict(case=b["case"],nearest_archived_case=int(ids[j]),distance=float(np.sqrt(max(0,1-sim[j]))),
                         closest_historical_leading_windows=win[ids[j]].tolist()))
    out=dict(historical_packet_sha256=SHA,archive_mode=mode,archive_parameters=len(tails),noncollapsed_clouds=len(hs),
             equal_points_per_cloud=tails.shape[1],historical_phase_screen=phases,comparisons=rows)
    write_json(results/"all-archive-comparison.json",out)
    print(len(tails),len(hs),len(phases),max([a["order"] for a in phases],default=0),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results",type=Path,default=RESULTS)
    p.add_argument("--archive",type=Path)
    a=p.parse_args()
    run(a.results,a.archive)

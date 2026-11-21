"""Replay the chaotic lift and independently extract the cubic root braid in two projections."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numba import njit
from model import apply_chart,charts,geometric_word,invariants,xyz


@njit(cache=True)
def itinerary(n,initial,D):
    t=initial;w=0j;theta=np.empty(n);u=np.empty(n,dtype=np.complex128)
    for i in range(n+768):
        w=.2*w+.6*D*np.exp(1j*t)
        t=(2*t+.12*np.sin(2*t)+.137)%(2*np.pi)
        if i>=768:theta[i-768]=t;u[i-768]=w
    return theta,u


def sample(row,folder,n):
    core=charts(row,max(4096,256*row["k"]))
    controls=[]
    for angle in (.229,.713):
        word=geometric_word(core,angle)
        linking,mu=invariants(word)
        assert linking==[0,0,0] and abs(mu)==row["abs_mu"],(row["case"],angle,linking,mu)
        controls.append(dict(angle=angle,geometric_crossings=len(word),pairwise_linking=linking,absolute_mu=abs(mu)))
    sets=[]
    for j,initial in enumerate((.27,1.34,3.16)):
        theta,u=itinerary(n,initial,row["return_residual_radius"])
        z=apply_chart(theta,u,row,j,core)
        sets.append(xyz(theta,z,row))
    file=folder/f"case-{row['case']}.npz"
    np.savez_compressed(file,points=np.concatenate(sets),component=np.repeat(np.arange(3,dtype=np.int8),n))
    return dict(case=row["case"],points=3*n,cloud_relative=f"frontiers/cubic/{file.name}",
                file_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),geometric_braid=controls)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificates",type=Path);ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--manifest",type=Path,required=True);ap.add_argument("--steps",type=int,default=28000)
    ap.add_argument("--cases",nargs="+",type=int);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=True);a.manifest.parent.mkdir(parents=True,exist_ok=True)
    rows=json.loads(a.certificates.read_text())["rows"]
    records=[sample(r,a.out,a.steps) for r in rows if a.cases is None or r["case"] in a.cases]
    a.manifest.write_text(json.dumps(records,indent=2)+"\n")
    print(f"{len(records)} clouds, two geometric braid checks per case.",flush=True)

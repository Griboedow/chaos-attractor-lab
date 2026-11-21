"""Measure cyclic order on independent tails and compare held-out phase labels."""
import argparse
import json
from pathlib import Path
from search import RESULTS, sample
from cyclic import scan_order
from phase_review import summarize
from validate import write_json
import numpy as np
from scipy.spatial import cKDTree


def oriented_gaps(z,m):
    groups=z[:len(z)-len(z)%m].reshape(-1,m,3).transpose(1,0,2)
    means=groups.mean(axis=1)
    centered=groups-means[:,None,:]
    _,bases=np.linalg.eigh(centered.transpose(0,2,1)@centered)
    coord=centered@bases
    low,high=coord.min(axis=1),coord.max(axis=1)
    half=(high-low)/2
    middle=means+np.einsum("mij,mj->mi",bases,(low+high)/2)
    all_gaps=[]
    for i in range(1,m):
        # A separating axis for oriented boxes is a face normal of either
        # box or a cross product of their edges. Unit axes retain native scale.
        aa=np.broadcast_to(bases[i].T,(i,3,3))
        bb=bases[:i].transpose(0,2,1)
        cross=np.cross(aa[:,:,None,:],bb[:,None,:,:]).reshape(i,9,3)
        directions=np.concatenate((aa,bb,cross),axis=1)
        norms=np.linalg.norm(directions,axis=2)
        directions/=np.maximum(norms,1e-14)[:,:,None]
        separation=np.abs(np.einsum("ijk,ik->ij",directions,middle[i]-middle[:i]))
        ai=np.abs(directions@bases[i])@half[i]
        aj=np.einsum("ijk,ik->ij",np.abs(directions@bases[:i]),half[:i])
        gaps=separation-ai-aj
        gaps[norms<1e-12]=-np.inf
        all_gaps.extend(gaps.max(axis=1).tolist())
    d=np.linalg.norm(np.ptp(z,axis=0))
    ar=np.array(all_gaps)/d
    return dict(finite_cloud_box_pairs=len(ar),pairs_with_positive_axis_gap=int((ar>1e-10).sum()),
                smallest_oriented_box_gap_fraction=float(ar.min()),median_oriented_box_gap_fraction=float(np.median(ar)))


def test(row,specified=None):
    z=np.load(row["cloud"])["points"].astype(float)
    if specified:
        m=specified
        v=summarize(z,m)["variance_explained"]
    else:
        m,v=scan_order(z[:36_000],48)
    if v<.55:
        return dict(case=row["case"],order=m,phase_variance_score=float(v),phase_separated=False)
    stats=summarize(z,m)
    orientation=oriented_gaps(z,m)
    cut=180_000
    tree=cKDTree(z[:cut])
    def classify(x):
        labels=tree.query(x,workers=1)[1]%m
        offsets=(labels-np.arange(len(x)))%m
        counts=np.bincount(offsets,minlength=m)
        return float(np.max(counts)/len(x))
    held=classify(z[cut:])
    starts=np.array([b["start"] for b in row["basin"][1:9]])
    w,tails=sample(np.tile(row["params"],(len(starts),1)),starts,120_000,180_000,6,24000)
    independent=[dict(start=x.tolist(),aligned_label_accuracy=classify(t),
                      rates=r.tolist()) for x,t,r in zip(starts,tails,w)]
    return dict(case=row["case"],phase_separated=True,**stats,**orientation,held_out_accuracy=held,
                independent_starts=independent)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results",type=Path,default=RESULTS)
    p.add_argument("--cases",type=int,nargs="*")
    p.add_argument("--orders",type=Path,help="Preview JSON specifying the tested higher-order partition")
    a=p.parse_args()
    orders={r["case"]:r["order"] for r in json.loads(a.orders.read_text())} if a.orders else {}
    for path in sorted((a.results/"checks").glob("*.json")):
        row=json.loads(path.read_text())
        if not row.get("cloud") or (a.cases and row["case"] not in a.cases):
            continue
        data=test(row,orders.get(row["case"]))
        write_json(a.results/"phase"/(path.name),data)
        print(row["case"],data["order"],data.get("minimum_box_gap_fraction"),data.get("held_out_accuracy"),flush=True)

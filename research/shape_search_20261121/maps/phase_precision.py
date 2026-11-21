"""Measure cyclic phase separation before binary32 export, across nine orbits."""
import argparse
import json
import sys
from pathlib import Path
from search import RESULTS, REPO, sample
from phase_validate import oriented_gaps
from validate import write_json
import numpy as np
from numba import njit
from scipy.spatial import cKDTree

sys.path.insert(0,str(REPO/"research"))
from twist_map import step_tangent


@njit(cache=True)
def independent_trace(p, start, burn, keep):
    x=start.copy()
    tangent=np.eye(3)
    z=np.empty((keep,3))
    for n in range(burn+keep):
        x,_=step_tangent(x,tangent,p)
        if n>=burn:
            z[n-burn]=x
    return z


def measure(results,case):
    row=json.loads((results/"checks"/f"{case}.json").read_text())
    phase=json.loads((results/"phase"/f"{case}.json").read_text())
    m=phase["order"]
    p=np.array(row["params"])
    start=np.array(row["primary_start"])
    _,native=sample(p[None,:],start[None,:],100_000,2_000_000,40,240_000)
    native=native[0]
    exported=np.load(row["cloud"])["points"]
    if not np.array_equal(native.astype("float32"),exported):
        raise RuntimeError("Native and retained clouds do not correspond")
    tree=cKDTree(native)
    def complete_periods(z,offset):
        begin=(-offset)%m
        n=(len(z)-begin)//m*m
        return z[begin:begin+n]
    alignments=[]
    kept=[complete_periods(native,0)]
    burn,keep=400_000,500_000
    for item in row["basin"][1:9]:
        start=np.array(item["start"])
        z=independent_trace(p,start,burn,keep)
        indices=np.arange(0,len(z),11)
        labels=tree.query(z[indices],workers=1)[1]%m
        counts=np.bincount((labels-indices)%m,minlength=m)
        offset=int(np.argmax(counts))
        accuracy=float(counts[offset]/len(indices))
        alignments.append(dict(start=start.tolist(),offset=offset,
            classified_points=len(indices),aligned_accuracy=accuracy))
        if accuracy>.99:
            kept.append(complete_periods(z,offset))
    combined=np.vstack(kept)
    write_json(results/"phase-binary64"/f"{case}.json",dict(case=case,phase_order=m,
        reference_binary64=oriented_gaps(native,m),export_rounding_max_coordinate_error=float(np.max(np.abs(native-exported))),
        independent_model="research/twist_map.py step_tangent",independent_burn=burn,
        independent_observed=keep,independent_starts=alignments,
        combined_binary64_points=len(combined),combined_binary64=oriented_gaps(combined,m)))
    print(case,m,len(combined),min(r["aligned_accuracy"] for r in alignments),flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results",type=Path,default=RESULTS)
    ap.add_argument("--cases",type=int,nargs="+",required=True)
    a=ap.parse_args()
    for case in a.cases:
        measure(a.results,case)

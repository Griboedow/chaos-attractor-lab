"""Plane-independent excursion detector using two fitted 3-D envelopes of the reference lobes."""
import argparse
import json
from pathlib import Path

import numpy as np
from numba import njit

from analyze import runs, write


def envelopes(z, rms):
    centers=[]; rotations=[]; low=[]; high=[]
    for p in (z[::2],z[1::2]):
        center=p.mean(0); axes=np.linalg.eigh(np.cov(p.T))[1]
        projected=(p-center)@axes
        centers.append(center); rotations.append(axes)
        low.append(projected.min(0)); high.append(projected.max(0))
    return tuple(np.asarray(a) for a in (centers,rotations,low,high))


@njit(cache=True)
def outside(z,centers,rotations,low,high,padding):
    missing=np.zeros(len(z),np.bool_)
    for i in range(len(z)):
        inside=False
        for b in range(2):
            dx=z[i,0]-centers[b,0];dy=z[i,1]-centers[b,1];dz=z[i,2]-centers[b,2]
            match=True
            for d in range(3):
                x=dx*rotations[b,0,d]+dy*rotations[b,1,d]+dz*rotations[b,2,d]
                if x<low[b,d]-padding or x>high[b,d]+padding:
                    match=False
                    break
            if match:
                inside=True
                break
        missing[i]=not inside
    return np.flatnonzero(missing)


def main(args):
    with np.load(args.cloud) as f:z=f["points"].astype(float)
    rms=np.sqrt(z.var(0).sum());model=envelopes(z[:125000],rms)
    cases=sorted(args.output.glob("cohort-*.json"))
    totals=[]
    for padding in (.05,.10):
        hold=outside(z[125000:],*model,padding*rms)
        cohorts=[]
        for case in cases:
            row=json.loads(case.read_text())
            trace=np.load(args.cache/row["trace"],mmap_mode="r")
            with np.load(args.cache/row["windows"]) as f: times=f["left_index"]
            ids=[]
            for a in range(0,len(trace),2000000):
                ids.append(a+outside(trace[a:a+2000000],*model,padding*rms))
            ids=np.concatenate(ids)
            parity_groups=runs(times,64); geometric_groups=runs(ids,64)
            p_to_g=sum(np.searchsorted(ids,b+18)>np.searchsorted(ids,a-16) for a,b,_ in parity_groups)
            g_to_p=sum(np.searchsorted(times,b+17)>np.searchsorted(times,a-16) for a,b,_ in geometric_groups)
            near=np.union1d(ids,ids-1)
            cohorts.append(dict(cohort=row["cohort"],observed=len(trace),outside_states=len(ids),geometric_groups=len(geometric_groups),
                         parity_groups=len(parity_groups), parity_events=len(times),parity_pairs_with_outside_endpoint=len(np.intersect1d(near,times)),
                         parity_groups_with_outside_state_near=int(p_to_g),geometric_groups_with_parity_pair_near=int(g_to_p)))
        totals.append(dict(padding_rms=padding,heldout_states=125000,heldout_outside=int(len(hold)),cohorts=cohorts))
    write(args.output/"envelopes.json",dict(description="Fit an oriented bounding box to each parity class in the first half of the native cloud, expand each face, then test every full-precision state geometrically against the union. Group gaps at <=64 and match groups with a sixteen-step margin.",model=dict(rms=rms,centers=model[0].tolist(),axes=model[1].tolist(),lower=model[2].tolist(),upper=model[3].tolist()),experiments=totals))


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    for name in ("cloud","cache","output"):
        p.add_argument("--"+name,type=Path,required=True)
    main(p.parse_args())

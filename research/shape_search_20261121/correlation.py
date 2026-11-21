"""Measure finite-scale pair correlations across independent starts of the five-component twist map."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/"maps"))
from search import sample
from index_campaign import RESULTS,write


def run(out):
    row=json.loads((RESULTS/"maps/checks/2503155.json").read_text())
    n=6
    rng=np.random.default_rng(2026112143)
    initial=rng.uniform(-.65,.65,(n,3))
    params=np.repeat([row["params"]],n,axis=0)
    rates,z=sample(params,initial,150000,800000,8,180000)
    assert rates.min()>.08
    rms=np.sqrt(np.mean(np.sum((z[0]-z[0].mean(0))**2,axis=1)))
    radii=np.geomspace(2e-7,1e-2,26)
    results=[]
    for i in range(0,n,2):
        left=z[i,rng.choice(z.shape[1],32000,replace=False)]
        right=z[i+1,rng.choice(z.shape[1],32000,replace=False)]
        a,b=cKDTree(left),cKDTree(right)
        counts=a.count_neighbors(b,radii*rms)
        corr=counts/(len(left)*len(right))
        mask=(radii>=1e-5)&(radii<=1e-3)
        assert np.min(counts[mask])>100
        slope=np.polyfit(np.log(radii[mask]),np.log(corr[mask]),1)[0]
        results.append({"first_start":i,"second_start":i+1,"cross_counts":counts.tolist(),
                        "cross_correlation":corr.tolist(),"slope_on_prescribed_band":float(slope)})
    write(out,{"case":2503155,"seed":2026112143,"initials":initial.tolist(),"native_precision":"binary64",
               "discarded":150000,"observed":800000,"sampled_per_start":32000,
               "native_rms":float(rms),"radii_in_rms_units":radii.tolist(),"fitted_radius_band":[1e-5,1e-3],
               "window_leading_exponents":rates.tolist(),"results":results,
               "scope":"Finite-scale correlation slope of sampled empirical measures across independent starts, distinct from the certified Hausdorff dimension of the attracting set."})
    print([r["slope_on_prescribed_band"] for r in results])


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,required=True)
    run(ap.parse_args().out)

"""Longer Lyapunov windows for cases with hundreds of visits per phase cycle."""
import argparse
import json
from search import RESULTS,START,sample
from validate import spectra,write_json
import numpy as np


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--cases",type=int,nargs="+",required=True)
    a=p.parse_args()
    for case in a.cases:
        row=json.loads((RESULTS/"checks"/f"{case}.json").read_text())
        pars=np.array(row["params"])[None,:]
        ref,_=sample(pars,START[None,:],500000,20_000_000,20,0)
        full,_,_,_,_=spectra(pars,START[None,:],300000,3_000_000,150000,0)
        starts=np.array([b["start"] for b in row["basin"]])
        rr,_=sample(np.tile(pars,(len(starts),1)),starts,150000,1_500_000,10,0)
        out=dict(case=case,reference_observed=20_000_000,reference_burn=500000,reference_windows=ref[0].tolist(),
                 independent_burn=300000,independent_observed=3_000_000,independent_windows=full[0].tolist(),
                 starts=starts.tolist(),cohort_burn=150000,cohort_observed=1_500_000,cohort_windows=rr.tolist(),
                 positive_all_long_windows=int((rr.min(axis=1)>0).sum()))
        write_json(RESULTS/"long-phase-controls"/f"{case}.json",out)
        print(case,float(ref.mean()),float(ref.min()),float(full[0,:,0].min()),out["positive_all_long_windows"],flush=True)

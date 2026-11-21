"""Regenerate a retained cloud directly from its saved parameters and start."""
import argparse
import hashlib
import json
from pathlib import Path
from search import CACHE,RESULTS,sample
import numpy as np


def regenerate(case,results,output):
    record=json.loads((results/"checks"/f"{case}.json").read_text())
    if not record.get("cloud"):
        raise ValueError("This case has no retained reference cloud")
    p=np.array(record["params"])[None,:]
    start=np.array(record["primary_start"])[None,:]
    rates,z=sample(p,start,100_000,2_000_000,40,240_000)
    pts=z[0].astype(np.float32)
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.resolve()==Path(record["cloud"]).resolve():
        raise ValueError("Choose a replay output separate from the archived cloud")
    np.savez_compressed(output,points=pts)
    sha=hashlib.sha256(output.read_bytes()).hexdigest()
    return dict(case=case,output=str(output.resolve()),points=len(pts),sha256=sha,
                archived_sha256=record["sha256"],byte_identical=(sha==record["sha256"]),
                leading_mean=float(rates.mean()),minimum_leading_window=float(rates.min()))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("case",type=int)
    p.add_argument("--results",type=Path,default=RESULTS);p.add_argument("--output",type=Path)
    a=p.parse_args();dest=a.output or CACHE/"replay"/f"{a.case}.npz"
    print(json.dumps(regenerate(a.case,a.results,dest),indent=2))

"""Check both expanding directions at multiple initial states on a longer run."""
import argparse
import json
from pathlib import Path
from search import RESULTS
from validate import spectra,write_json
import numpy as np


def verify(case,results):
    row=json.loads((results/"checks"/f"{case}.json").read_text())
    starts=np.array([b["start"] for b in row["basin"]])
    p=np.tile(row["params"],(len(starts),1))
    w,_,_,_,_=spectra(p,starts,100_000,500_000,25_000,0)
    rows=[]
    for start,windows in zip(starts,w):
        rows.append(dict(start=start.tolist(),mean=windows.mean(axis=0).tolist(),
            min_window=windows.min(axis=0).tolist(),max_window=windows.max(axis=0).tolist()))
    packet=dict(case=case,burn=100_000,observed=500_000,window=25_000,starts=rows,
                all_second_windows_positive=sum(a["min_window"][1]>0 for a in rows),
                maximum_volume_identity_error=float(np.max(abs(w.sum(axis=2)-3*np.log(p[0,0])))))
    write_json(results/"two-rate-controls"/f"{case}.json",packet)
    print(case,packet["all_second_windows_positive"],min(r["min_window"][1] for r in rows),flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results",type=Path,default=RESULTS)
    a=ap.parse_args()
    for case in json.loads((a.results/"proposals-hyper.json").read_text()):
        verify(int(case),a.results)

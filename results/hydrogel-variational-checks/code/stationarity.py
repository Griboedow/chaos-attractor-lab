"""Fixed-window tip statistics; a drift diagnostic, not a stationarity test."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('folders',nargs='+',type=Path)
ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();out=[]
for folder in args.folders:
    with np.load(folder/'trajectory.npz') as d:t=d['t'];tip=d['tip']
    rows=[]
    for left in np.arange(t[0],t[-1]-1e-8,24):
        right=min(left+24,t[-1]);ids=(t>=left)&(t<right)
        block=tip[ids]
        rows.append(dict(start=float(left),end=float(right),samples=len(block),
            mean=block.mean(axis=0).tolist(),std=block.std(axis=0).tolist(),
            minimum=block.min(axis=0).tolist(),maximum=block.max(axis=0).tolist()))
    row=dict(run=folder.name,windows=rows,
        input_sha256=hashlib.sha256((folder/'trajectory.npz').read_bytes()).hexdigest())
    out.append(row)
    print(json.dumps(dict(run=folder.name,first_height_mean=rows[0]['mean'][2],
        last_height_mean=rows[-1]['mean'][2],last_height_std=rows[-1]['std'][2])))
args.output.write_text(json.dumps(dict(window_duration=24,runs=out,
    scope='Unaligned physical tip coordinates over consecutive windows; descriptive drift diagnostics only.'),indent=2)+'\n',encoding='utf-8')

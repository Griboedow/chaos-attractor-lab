"""Matched phase seeds from the previous chaotic trajectory and from global initial states."""
import argparse,json
import numpy as np
from native_section import CASES,OLD,RAW,orbit_section,raw_trace

def starts():
    ref=np.load(CASES['coexisting_chaotic_cloud']['path'])['points']
    rng=np.random.default_rng(564720261121)
    observed=rng.choice(len(ref),64,replace=False)
    rs=[dict(kind='cloud',source_index=int(i),initial=ref[i].tolist()) for i in observed]
    basins=json.loads((OLD/'basins.json').read_text(encoding='utf-8'))['broad']
    rs += [dict(kind='broad',source_index=i,initial=v['initial']) for i,v in enumerate(basins)]
    return rs

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--step',type=float,required=True);ap.add_argument('--duration',type=float,default=60000)
    ap.add_argument('--indices',type=int,nargs='+',required=True);a=ap.parse_args()
    inputs=starts()
    for i in a.indices:
        v=inputs[i]
        path=RAW/f"s{i:03d}-h{a.step:g}-t{a.duration:g}.npz"
        raw_trace(path,v['initial'],a.step,a.duration)

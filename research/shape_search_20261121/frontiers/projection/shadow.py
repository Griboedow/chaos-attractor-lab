"""Restart independent MPFR dynamics at stored pre-event states to audit local itinerary."""
import argparse
import json
from pathlib import Path
import sys

import gmpy2 as g
import numpy as np

EXTENSION = Path(__file__).resolve().parents[2] / "extension" / "projection"
sys.path.insert(0, str(EXTENSION))
from precision import metric_step


def main(args):
    doc = json.loads(args.record.read_text())
    planes = json.loads(args.planes.read_text())
    n = np.asarray(planes["normals"]["frozen"])
    cut = planes["cuts"]["frozen"]
    with np.load(args.events) as f:
        windows = f["xyz"]
        indices = f["left_index"]
        offsets = f["offsets"]
    out = []
    with g.context(precision=args.bits):
        step = metric_step(doc["parameters"])
        mn = tuple(map(g.mpfr, n))
        mc = g.mpfr(cut)
        for k, segment in zip(indices, windows):
            state = tuple(map(g.mpfr, segment[0]))
            sides = []
            errs = []
            for j in range(len(segment)):
                if j:
                    state = step(state)
                sides.append(sum(u*v for u,v in zip(mn,state)) > mc)
                errs.append(np.linalg.norm(np.asarray([float(v) for v in state]) - segment[j]))
            native = segment @ n > cut
            left = int(np.flatnonzero(offsets == 0)[0])
            out.append((int(k), bool(sides[left] == sides[left+1]), np.count_nonzero(np.asarray(sides) != native),
                        float(max(errs)), float(max(errs[left:left+2]))))
    rows = np.array(out)
    args.out.write_text(json.dumps(dict(bits=args.bits, case=doc["case"], initial_offset=int(offsets[0]), final_offset=int(offsets[-1]),
        windows=len(rows), event_reproduced=int(rows[:,1].sum()), full_itinerary_matches=int((rows[:,2] == 0).sum()),
        max_endpoint_euclidean_error=float(rows[:,4].max()), max_window_euclidean_error=float(rows[:,3].max()),
        endpoint_error_quantiles=np.quantile(rows[:,4], [0,.25,.5,.75,1]).tolist(), mismatching_events=rows[rows[:,1]==0,0].astype(int).tolist(),
        method="Restart each MPFR trajectory from the exact binary64 state sixteen steps before its event. Compare the independently propagated signs up to offset +17."),indent=2)+"\n",encoding="utf8")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ("record","planes","events","out"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--bits",type=int,default=256)
    main(p.parse_args())

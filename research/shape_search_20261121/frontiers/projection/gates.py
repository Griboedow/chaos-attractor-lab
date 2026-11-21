"""Compare handoff points from the dense reference lobes to matched random states."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from analyze import geometry, median_range, write


def main(args):
    with np.load(args.cloud) as f: source=f["points"].astype(float)
    docs=json.loads((args.output/"episodes.json").read_text())
    rec=json.loads(args.record.read_text())
    planes=json.loads((args.output/"planes.json").read_text())
    pl=np.array(planes["normals"]["frozen"]); cut=planes["cuts"]["frozen"]
    rms=planes["rms"]; limit=docs["cloud_distance_threshold_rms"]
    tree=cKDTree(source[:125000]); rng=np.random.default_rng(229866505)
    extracts=[]
    for s in sorted(args.output.glob("cohort-*.json")):
        stat=json.loads(s.read_text());cohort=stat["cohort"]
        traj=np.load(args.cache/stat["trace"],mmap_mode="r")
        for event in (p for p in docs["episodes"] if p["cohort"]==cohort):
            a=event["first"]; b=event["last"]
            x=traj[a-4096:a+1]
            dist=tree.query(x)[0]/rms
            quiet=dist<=limit
            intact=np.convolve(quiet.astype(int),np.ones(32,dtype=int),"valid")==32
            j=int(np.flatnonzero(intact)[-1]+31)
            assert j<len(x)-1
            pre=a-4096+j
            enter=traj[pre];out=traj[pre+1]; path=traj[pre-8:pre+9]
            yz=geometry(path,rec["parameters"])
            # Last quiet state before the first pair, then the prediction for one step.
            extracts.append(dict(cohort=cohort,episode_start=a,episode_end=b,quiet_end=pre,
                steps_to_first_pair=a-pre,xyz=enter.tolist(),next_xyz=out.tolist(),side=int(enter@pl>cut),
                geometry=yz[8].tolist(),next_geometry=yz[9].tolist(),distance_rms=dist[j],next_distance_rms=dist[j+1]))
    discovery=[v for v in extracts if v["cohort"]<=16]
    blind=[v for v in extracts if v["cohort"]>16]
    # Calibrate the nearest-discovery-gate score against 10k random states
    # conditionally matched to the side of every individual held-out departure.
    for t in blind:
        search=[v for v in discovery if v["side"]==t["side"]]
        if not search:
            continue
        dn=cKDTree(np.array([v["xyz"] for v in search]))
        held=source[125000:]
        pool=held[(held@pl>cut)==t["side"]]
        baseline=dn.query(pool[rng.integers(0,len(pool),10000)])[0]
        distance=dn.query(t["xyz"])[0]
        t["distance_to_discovery_gate_rms"]=distance/rms
        t["fraction_of_matched_reference_states_as_close"]=(1+int((baseline<=distance).sum()))/(1+len(baseline))
    scores=[p["fraction_of_matched_reference_states_as_close"] for p in blind if "fraction_of_matched_reference_states_as_close" in p]
    write(args.output/"gates.json",dict(criterion="Last point of the last 32 consecutive states at or below the held-out 99.9% dense-cloud distance cut, within 4096 steps before the first violation of each 64-gap episode.",
       discovery_episodes=len(discovery),new_initial_condition_episodes=len(blind),heldout_percentile_score_quantiles=median_range(scores),
       median_handoff_lead_steps=median_range([v["steps_to_first_pair"]for v in extracts]),rows=extracts))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("record","cloud","cache","output"):
        p.add_argument("--"+name,required=True,type=Path)
    main(p.parse_args())

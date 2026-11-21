"""Test candidate gateway neighborhoods against every state of eight held-out traces."""
import argparse
import json
from pathlib import Path

import numpy as np
from numba import njit
from scipy.spatial import cKDTree

from analyze import runs, write


@njit(cache=True)
def gate_distances(z,gates,low,high,maxdist):
    values=np.full(len(z),np.inf)
    for i in range(len(z)):
        p=z[i]
        if p[0]<low[0] or p[0]>high[0] or p[1]<low[1] or p[1]>high[1] or p[2]<low[2] or p[2]>high[2]:
            continue
        best=maxdist*maxdist
        for g in gates:
            dx=p[0]-g[0];dy=p[1]-g[1];dz=p[2]-g[2]
            d=dx*dx+dy*dy+dz*dz
            if d<best:best=d
        if best<maxdist*maxdist:values[i]=np.sqrt(best)
    return values


def main(args):
    planes=json.loads((args.output/"planes.json").read_text())
    gates=json.loads((args.output/"gates.json").read_text())["rows"]
    episodes=json.loads((args.output/"episodes.json").read_text())["episodes"]
    limits=[.01,.03,.06,.10]
    known=np.array([v["xyz"] for v in gates if v["cohort"]<=16])
    rms=planes["rms"]; reach=max(limits)*rms; low=known.min(0)-reach; high=known.max(0)+reach
    with np.load(args.cloud) as f:reference=f["points"].astype(float)
    tree=cKDTree(reference[:125000]);gate_tree=cKDTree(known)
    cf=json.loads((args.output/"episodes.json").read_text())["cloud_distance_threshold_rms"]
    normal=np.array(planes["normals"]["frozen"]);cut=planes["cuts"]["frozen"]
    audit=[]; controls=[]; rng=np.random.default_rng(9058739)
    for cohort in range(17,25):
        record=json.loads((args.output/f"cohort-{cohort:02d}.json").read_text())
        z=np.load(args.cache/record["trace"],mmap_mode="r")
        bursts=[p for p in episodes if p["cohort"]==cohort]
        first=np.array([p["first"] for p in bursts])
        indexes=[[] for _ in limits]
        for i in range(0,len(z),2000000):
            dist=gate_distances(z[i:i+2000000],known,low,high,reach)/rms
            for row,d in zip(indexes,limits):row.append(i+np.flatnonzero(dist<d))
        for d,row in zip(limits,indexes):
            times=np.concatenate(row); alarms=runs(times,32)
            # Prospective alarm time is the FIRST hit. It is never moved later
            # using the subsequently observed parity event.
            alarm=np.array([a for a,b,n in alarms],dtype=int)
            predicted=sum(np.searchsorted(first,a+201)>np.searchsorted(first,a+1) for a in alarm)
            captured=sum(np.searchsorted(alarm,b)>np.searchsorted(alarm,b-200) for b in first)
            audit.append(dict(cohort=cohort,radius_rms=d,foreground_states=len(times),alarms=len(alarm),
                      alarms_followed_by_episode_within_200=int(predicted),episode_starts=len(first),episodes_with_prior_alarm=int(captured)))
        # Compare ordinary density-boundary crossings away from any episode
        # with actual departures, applying the identical 32-in / one-out rule.
        ids=rng.integers(0,len(z)-4096,192)
        for i in ids:
            if any(p["last"]+4096>=i and p["first"]-4096<=i+4096 for p in bursts):continue
            x=z[i:i+4096];dist=tree.query(x)[0]/rms
            ok=dist<=cf
            valid=(np.convolve(ok.astype(int),np.ones(32,dtype=int),"valid")[:-1]==32)&~ok[32:]
            xs=x[np.flatnonzero(valid)+31]
            if len(xs):
                ds=gate_tree.query(xs)[0]/rms
                controls.extend((cohort,bool(xx@normal>cut),float(dd)) for xx,dd in zip(xs,ds))
    counts=[]
    for d in limits:
        rows=[v for v in audit if v["radius_rms"]==d]
        counts.append(dict(radius_rms=d,**{k:sum(v[k]for v in rows)for k in ("foreground_states","alarms","alarms_followed_by_episode_within_200","episode_starts","episodes_with_prior_alarm")}))
    depart=np.array([p["distance_to_discovery_gate_rms"] for p in gates if p["cohort"]>=17])
    cd=np.array([v[2] for v in controls if not v[1]])
    comparison=[dict(radius_rms=d,true_departures_inside=int((depart<d).sum()),ordinary_negative_lobe_crossings_inside=int((cd<d).sum()))for d in limits]
    write(args.output/"prospective.json",dict(method="Learn 17 departure-point locations from original and prior cohorts, scan all 200 million later-cohort states for first entry into balls, merge gaps <=32, predict an episode onset in the next 1..200 steps. Radius comparisons are exploratory.",
        totals=counts,by_cohort=audit,ordinary_crossing_controls=dict(sampled_blocks_per_cohort=192,block_length=4096,
        exclusion_margin=4096,sampled_boundary_crossings=len(controls),negative_side_crossings=len(cd),positive_side_crossings=len(controls)-len(cd),
        true_departures=len(depart),comparison=comparison)))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("cache","cloud","output"):
        p.add_argument("--"+name,required=True,type=Path)
    main(p.parse_args())

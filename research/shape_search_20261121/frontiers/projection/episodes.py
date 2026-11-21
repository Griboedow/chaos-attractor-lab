"""Resolve rare excursions in state space and test whether their timing depends on the plane."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from analyze import event_times, fit, geometry, make_planes, median_range, runs, sha, write


def intervals(ids, expansion):
    return [(max(0,a-expansion),b+expansion) for a,b,_ in runs(ids,64)]


def covered(groups, times):
    return sum(any(a <= t <= b for t in times[np.searchsorted(times,a):np.searchsorted(times,b+1)]) for a,b in groups)


def main(args):
    record=json.loads(args.record.read_text())
    ref=json.loads(args.reference.read_text())
    with np.load(args.cloud) as f: source=f["points"].astype(float)
    tree=cKDTree(source)
    planes=make_planes(source,ref); normal,cut=planes[0][1:]
    rms=ref["source_rms"]
    steps=np.array([-4096,-1024,-256,-64,-32,-16,-8,-4,-2,-1,0,1,2,4,8,16,32,64,256,1024,4096])
    summaries=[]; phase=[]; analysis=[]
    train,hold=source[:125000],source[125000:]
    trained=cKDTree(train)
    # Fixed geometric out-of-support score, calibrated on a separate reference half.
    held_dist=trained.query(hold)[0]/rms
    threshold=float(np.quantile(held_dist,.999))
    gate_model=None
    for summary in sorted(args.output.glob("cohort-*.json")):
        v=json.loads(summary.read_text())
        path=args.cache/v["trace"]
        z=np.load(path,mmap_mode="r")
        with np.load(args.cache/v["windows"]) as f: times=f["left_index"]
        clusters=runs(times,64)
        if not clusters: continue
        starts=np.array([item[0] for item in clusters])
        ends=np.array([item[1] for item in clusters])
        assert starts.min()>4096 and ends.max()+4096<len(z)
        onset=z[starts[:,None]+steps]
        recovery=z[ends[:,None]+steps]
        rms_on=trained.query(onset.reshape(-1,3))[0].reshape(len(starts),len(steps))/rms
        rms_off=trained.query(recovery.reshape(-1,3))[0].reshape(len(ends),len(steps))/rms
        focus = 256
        for j,(a,b,n) in enumerate(clusters):
            w = z[a-focus:b+focus+1]
            score=trained.query(w)[0]/rms
            target=w[focus:focus+b-a+2]
            near=geometry(target,record["parameters"])
            flag=score>threshold
            inside=flag[focus:focus+b-a+2]
            before=np.flatnonzero(~flag[:focus]); after=np.flatnonzero(~flag[focus+b-a+2:])
            scope=scope_last=None
            if flag[focus] and len(before):scope=a-focus+int(before[-1])+1
            if flag[focus+b-a+1] and len(after):scope_last=b+1+int(after[0])
            on_g=geometry(onset[j],record["parameters"])
            xy0=z[a]; xy1=z[b+1]
            summaries.append(dict(cohort=v["cohort"],first=a,last=b,total=n,span=b-a+1,phase_flipped=bool(n%2),
                    positive_count=int(np.sum(z[times[(times>=a)&(times<=b)]]@normal>cut)),
                    start_xyz=xy0.tolist(),end_xyz=xy1.tolist(),
                    nearest_rms_on=rms_on[j].tolist(),nearest_rms_off=rms_off[j].tolist(),
                    off_support_fraction=float(inside.mean()), maximum_cloud_dist_rms=float(score.max()),
                    min_axis_in_episode=float(near[:,0].min()),min_core_in_episode=float(near[:,1].min()),
                    onset_core_distance=on_g[:,1].tolist(), before_first_offsupport=scope,after_last_offsupport=scope_last))
        if v["cohort"]==0:
            nominal=intervals(times,16)
            for label,n,cutoff in planes:
                other=times if label=="frozen" else event_times(z,n,cutoff)
                other_intervals=intervals(other,16)
                analysis.append(dict(plane=label,episodes=len(other_intervals),
                        original_episode_matched=covered(nominal,other),alt_episode_matched=covered(other_intervals,times)))
        phase.extend((v["cohort"],a,b,n) for a,b,n in clusters)
    primary=[e for e in summaries if e["cohort"]==0]
    others=[e for e in summaries if e["cohort"]!=0]
    def group(rows):
        if not rows:return {}
        return dict(episodes=len(rows),distinct_starts=len(set(q["cohort"]for q in rows)),flipped=sum(q["phase_flipped"]for q in rows),
               median_nearest_rms_on=np.median([q["nearest_rms_on"] for q in rows],axis=0).tolist(),
               median_nearest_rms_off=np.median([q["nearest_rms_off"] for q in rows],axis=0).tolist(),
               median_onset_core_distance=np.median([q["onset_core_distance"]for q in rows],axis=0).tolist(),
               off_support_fraction=median_range([q["off_support_fraction"] for q in rows]),
               maximum_distance_rms=median_range([q["maximum_cloud_dist_rms"] for q in rows]))
    write(args.output/"episodes.json",dict(coordinates="native XYZ; no centering or PCA in distances",time_offsets=steps.tolist(),
        calibration="Tree on first 125000 reference states; cutoff at 99.9th percentile of nearest distances of second 125000.",
        cloud_distance_threshold_rms=threshold,primary=group(primary),independent=group(others),
        primary_plane_sensitivity=analysis,episodes=summaries))
    print(len(primary),len(others),threshold,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for field in ("record","reference","cloud","cache","output"):
        p.add_argument("--"+field,type=Path,required=True)
    main(p.parse_args())

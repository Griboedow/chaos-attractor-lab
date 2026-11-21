"""Verify all windows against complete traces and produce compact portable audit indexes."""
import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from analyze import sha, write


def main(args):
    docs=[json.loads(p.read_text()) for p in sorted(args.output.glob("cohort-*.json"))]
    assert [r["cohort"] for r in docs]==list(range(25))
    eps=json.loads((args.output/"episodes.json").read_text())["episodes"]
    center=json.loads((args.output/"planes.json").read_text())
    index=[]; allids={}
    for v in docs:
        actual=args.cache/v["trace"]; samples=args.cache/v["windows"]
        z=np.load(actual,mmap_mode="r")
        with np.load(samples) as w:
            t=w["left_index"]
            assert np.array_equal(z[t[:,None]+w["offsets"]],w["xyz"])
            allids[str(v["cohort"])]=t.tolist()
        index.append(dict(cohort=v["cohort"],states=v["observed"],events=v["events"],windows_sha256=sha(samples),
             windows_relative="frontiers/projection/"+samples.name,trace_relative="frontiers/projection/"+actual.name,
             record=f"cohort-{v['cohort']:02d}.json"))
    with ThreadPoolExecutor(max_workers=2) as pool:
        for row,checksum in zip(index,pool.map(sha,[args.cache/v["trace"]for v in docs])):
            row["trace_sha256"]=checksum
    write(args.output/"manifest.json",dict(index_base="shape-search-20261121 cache root",trace_format="float64 N by 3 native XYZ in zero-based post-burn order",
        window_format="float64 native XYZ; left_index + offsets -16 through +17 inclusive; both endpoints included",
        source_cloud_sha256=center["cloud_sha256"],records=index))
    write(args.output/"indices.json",dict(description="Zero-based left endpoint k of an equal-label pair at states k,k+1 after 100000 burn-in iterations, beginning with x_100001.",
         cohorts=allids))
    counts=[]
    for label,ids in [("original",[0]),("earlier_other",list(range(1,17))),("new_starts",list(range(17,25))),("all",list(range(25)))]:
        subset=[v for v in docs if v["cohort"] in ids]; seq=[v for v in eps if v["cohort"] in ids]
        counts.append(dict(cohort=label,starts=len(subset),represented_starts=sum(bool(v["events"]) for v in subset),
               states=sum(v["observed"]for v in subset),pairs_compared=sum(v["compared"]for v in subset),
               events=sum(v["events"]for v in subset),episodes=len(seq),phase_changes=sum(v["phase_flipped"]for v in seq),
               mean_episode_count_per_million=len(seq)/sum(v["observed"]for v in subset)*1e6))
    write(args.output/"summary.json",dict(case=2298665,groups=counts,unique_event_windows=sum(v["events"]for v in docs),
        all_saved_windows_equal_complete_native_trace=True,manifest_sha256=sha(args.output/"manifest.json"),
        indices_sha256=sha(args.output/"indices.json")))
    print(counts,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("cache","output"):
        p.add_argument("--"+name,required=True,type=Path)
    main(p.parse_args())

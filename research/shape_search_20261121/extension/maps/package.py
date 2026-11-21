"""Audit coefficient duplicates, reproduce clouds, and package replayable certificates."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
from multistable import ROOT,OLD
import numpy as np
from cycle_tubes import verify
from search import START,sample,design
from phase_refine import local_design,BASES
from validate import spectra

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def get(path):
    return json.loads(path.read_text())

def pack(cache,out):
    oc=cache.parents[1]/"maps"
    inputs=[np.load(f)["params"] for pattern in ("batch-*/scout.npz","hyper/batch-*/scout.npz") for f in sorted(oc.glob(pattern))]
    for f in sorted(oc.glob("cycle-*/screen.json")):
        b=get(f);inputs.append(design(b["batch"],b["count"]))
    for pattern in ("nested-phase/screen-*.json","nested-deep/screen-*.json"):
        for f in sorted(oc.glob(pattern)):
            b=get(f);g=next(i for i,(case,_) in enumerate(BASES) if case==b["base"])
            p=np.array(get(OLD/"checks"/(str(b["base"])+".json"))["params"])
            inputs.append(local_design(p,g,b["parameters"]))
    old={tuple(row) for x in inputs for row in x}
    new=[tuple(row) for f in sorted(cache.glob("design-*.npz")) for row in np.load(f)["params"]]
    if len(set(new))!=1920 or old.intersection(new):raise ValueError("Duplicate design")
    results=[];old_final={r["case"] for r in get(OLD/"manifest.json")["finalists"]}
    (cache/"clouds").mkdir(exist_ok=True)
    for f in sorted((out/"traps").glob("*.json")):
        item=get(f);box=item["box"];case=box["case"];m=box["case_order"];p=np.array(box["params"])
        cloud=cache/"clouds"/f"{case}.npz"
        if not cloud.exists():
            w,z=sample(p[None,:],START[None,:],100000,2000000,40,240000)
            np.savez_compressed(cloud,points=z[0].astype(np.float32))
        q=np.load(cloud)["points"]
        if q.shape!=(240000,3) or not np.isfinite(q).all():raise ValueError("Bad cloud")
        ck=out/"checks"/f"{case}.json"
        if ck.exists():mle=get(ck)["independent_spectrum"][0]
        else:
            windows,*_=spectra(p[None,:],START[None,:],80000,600000,30000,32768)
            mle=float(windows[0,:,0].mean())
        nc=out/"noncollapse"/f.name
        result=get(nc) if nc.exists() else {}
        results.append(dict(case=case,family="Dissipative three-axis twist map",label=f"{m} cyclic continua",components=m,
            previously_selected=case in old_final,footprint="reanalysis of previously screened coefficients",key="points",
            cloud=str(cloud.resolve()),cloud_relative=cloud.relative_to(cache.parents[1]).as_posix(),sha256=digest(cloud),
            independent_finite_time_largest_exponent=mle,noncollapse=result.get("distinct_second_return_orbit"),
            canonical_start_entry_step=result.get("canonical_start_entry_step"),trap_file="traps/"+f.name,
            proof_sha256=digest(f),points=240000))
    cycle_dir=out/"cycle-certificates";cycle_dir.mkdir(exist_ok=True)
    controls=[]
    for hit in get(out/"cycles.json"):
        source=cache/hit["certificate"];p=get(source);dest=cycle_dir/source.name
        if not all(verify(p["params"],c) for c in p["cycles"]):raise ValueError("Cycle replay failed "+p["id"])
        shutil.copyfile(source,dest)
        controls.append(dict(**hit,sha256=digest(dest),relative_certificate="cycle-certificates/"+dest.name))
    summary=dict(archived_unique_designs_checked=len(old),new_unique_designs=len(set(new)),overlap_with_archive=len(old.intersection(new)),
        archived_high_phase_candidates_rechecked=24,new_designs_with_scanned_stable_period_eight_or_larger=len(get(cache/"scout.json")),
        new_vectors_with_certified_sinks=len(controls),new_vectors_with_proven_coexistence=sum(len(r["periods"])>1 for r in controls),
        finalist_count=len(results),finalists=results,cycles=controls)
    (out/"manifest.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("Certified components",[(r["case"],r["components"],r["noncollapse"]) for r in results],flush=True)
    print("Cycle packets",len(controls),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("cache",type=Path);ap.add_argument("out",type=Path)
    a=ap.parse_args();pack(a.cache,a.out)

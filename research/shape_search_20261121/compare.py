"""Compare saved 3-D samples with the frozen catalog and archived README states."""
import argparse
import gzip
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "research"))
sys.path.insert(0, str(ROOT / "tools"))
from compare_affine_shapes import whiten
from compare_attractor_shapes import align, distance, fingerprint, normalize
from data_cache import data_path, records
from render_portraits import coordinates
from animate_portraits import load_views


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(obj,indent=2,allow_nan=False)+"\n").encode("utf-8")
    path.write_bytes(gzip.compress(payload,mtime=0) if path.suffix==".gz" else payload)


def downsample(z, n=24000, seed=17):
    assert z.ndim==2 and z.shape[1]==3,z.shape
    if len(z) > n:
        idx = np.random.default_rng(seed).choice(len(z),n,replace=False)
        z = z[np.sort(idx)]
    return np.asarray(z,dtype=np.float64)


def digest_file(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def resolve_source(name):
    if name.startswith("external/"):
        return Path(os.environ.get("CHAOS_SOURCE_ARCHIVE","D:/Codex")) / name.removeprefix("external/")
    return data_path(name)


def original_portrait(family,case):
    if family=="fastica":
        sub=f"portraits-01/case-{case}.npz"
        records= json.loads((ROOT/"results/fastica/external-data.json").read_text(encoding="utf-8"))
        sha=next(r["sha256"] for r in records if r["path"]==sub)
        name="external/unexpected-attractors-20261101/fastica/"+sub
        path=resolve_source(name)
        assert digest_file(path)==sha,name
        with np.load(path) as z:
            # These are the three physical columns of the published portrait.
            return name,downsample(z["cloud"][:,[0,1,3]])
    folder="yogi-geometry" if case=="4214" else "yogi-plane"
    inventory=json.loads((ROOT/"results"/folder/"external-inventory.json").read_text(encoding="utf-8"))
    receipt=next(r for r in inventory["files"] if Path(r["path"].replace("\\","/")).name==f"case-{case}.npy")
    name="external/"+receipt["path"].replace("\\","/").split("Codex/",1)[1]
    path=resolve_source(name)
    assert digest_file(path)==receipt["sha256"],name
    return name,downsample(np.load(path,mmap_mode="r"))


def reference_inputs():
    prefix = "results/attractor-novelty/reference/trajectories/"
    paths = {p.relative_to(ROOT).as_posix() for p in (ROOT/prefix).glob("*.npz")}
    paths.update(k for k in records() if k.startswith(prefix))
    assert len(paths) == 108
    for name in sorted(paths):
        with np.load(data_path(name)) as f:
            yield "reference/"+Path(name).stem,name,downsample(f["points"])
    catalog = json.loads((ROOT/"records/portrait-catalog.json").read_text(encoding="utf-8"))
    for row in catalog["rows"]:
        if row["id"].startswith("shapes-20261121-"):
            continue
        if row.get("featured") and all(p["mode"]=="histogram" for p in row["panels"]):
            views,receipts = load_views(row)
            for panel,points,receipt in zip(row["panels"],views,receipts):
                yield row["id"]+"/"+panel["label"],receipt["source"],downsample(points)
            continue
        for panel in row["panels"]:
            if row["id"] in ("fastica","yogi-expansion"):
                source,z=original_portrait(row["id"],panel["label"])
                yield row["id"]+"/"+panel["label"],source,z
            elif panel["mode"] == "points":
                ident = row["id"]+"/"+panel["label"]
                yield ident,panel["source"],downsample(coordinates(panel))
            elif row["id"] in ("tetrahedral","raar-constellation-columns","raar-constellation-fans",
                               "rk23","adamax","yogi","yogi-arc"):
                ident = row["id"]+"/"+panel["label"]
                number=panel["label"]
                folder=Path(panel["source"]).parent
                key="states"
                if row["id"]=="tetrahedral":
                    source=f"results/unexpected-search/tetra-render/sample-{number}.npz"
                    key="points"
                elif row["id"].startswith("raar-constellation-"):
                    source=(folder/f"case-{number}-states.npz").as_posix()
                elif row["id"]=="yogi-arc":
                    source=(folder/"animation.npz").as_posix()
                else:
                    source=(folder/f"case-{number}-animation.npz").as_posix()
                with np.load(data_path(source)) as f:
                    yield ident,source,downsample(f[key])
    for family,name in (("archive/twist","results/twist-geometry/long-runs.npz"),
                        ("archive/vortex","results/gaussian-vortex/long-runs.npz")):
        with np.load(data_path(name)) as f:
            for index,(case,z) in enumerate(zip(f["ids"],f["tail"])):
                yield f"{family}/{case}/{index}",name,downsample(z)


def descriptors(z):
    outputs = {}
    for mode,transform in (("similarity",normalize),("affine",whiten)):
        try:
            shape = transform(z)
        except ValueError:
            continue
        outputs[mode] = fingerprint(shape)
    return outputs


def split_sample(key,z):
    if key.split("/")[0] in ("cycles","links","skeletons"):
        # A closed polygon can represent one traversal only. Its two time halves
        # are different arcs, so distribute both control samples over the orbit.
        return z[::2],z[1::2],"Interleaved vertices of the closed sampled geometry"
    return z[:len(z)//2],z[len(z)//2:],"Disjoint temporal halves"


def prepare_refs(out):
    rows, hashes = [], {}
    for name, source, z in reference_inputs():
        if source not in hashes:
            hashes[source] = digest_file(resolve_source(source))
        fps = descriptors(z)
        rows.append({"name":name,"source":source,"sha256":hashes[source],"sample_count":len(z),
                     "sample_sha256":hashlib.sha256(z.astype("<f8").tobytes()).hexdigest(),
                     "descriptors":{k:v.tolist() for k,v in fps.items()}})
    write(out,{"reference_count":len(rows),"rows":rows})
    print(f"{len(rows)} labeled reference samples",flush=True)


def compare(index, ref_file, cache, dest, fit):
    references = json.loads(ref_file.read_text(encoding="utf-8"))["rows"]
    items = json.loads(index.read_text(encoding="utf-8"))
    candidates = {}
    for key,source in items.items():
        with np.load(cache/source) as data:
            z = downsample(data["points"])
        candidates[key] = (z,descriptors(z))
    result, needed = {},set()
    for key,(z,fps) in candidates.items():
        half1,half2,protocol=split_sample(key,z)
        row = {"cloud":items[key], "compared":{},"split_protocol":protocol}
        for mode,fp in fps.items():
            transform = normalize if mode == "similarity" else whiten
            a,b = transform(half1),transform(half2)
            selfdistance = float(distance(fingerprint(a),fingerprint(b,1)))
            near = sorted((float(distance(fp,ref["descriptors"][mode])),ref["name"])
                          for ref in references if mode in ref["descriptors"])
            cohort = sorted((float(distance(fp,other_fps[mode])),other)
                            for other,(_,other_fps) in candidates.items()
                            if other != key and mode in other_fps)
            row["compared"][mode] = {"split_distance":selfdistance,"reference_count":len(near),
                                     "nearest_reference":near[:3],"nearest_candidate":cohort[:3]}
            needed.add(near[0][1])
        result[key] = row
    if fit:
        raw = {n:z for n,_,z in reference_inputs() if n in needed}
        for key,row in result.items():
            z = candidates[key][0]
            for mode,test in row["compared"].items():
                transform = normalize if mode == "similarity" else whiten
                nearest = test["nearest_reference"][0][1]
                a,b = downsample(transform(z),1000),downsample(transform(raw[nearest]),1000,19)
                half1,half2,_=split_sample(key,z)
                selfa = downsample(transform(half1),1000)
                selfb = downsample(transform(half2),1000,19)
                test["local_icp_to_nearest"],_ = align(a,b)
                test["local_icp_between_halves"],_ = align(selfa,selfb)
            write(dest, result)
    write(dest,result)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference",type=Path,required=True)
    ap.add_argument("--prepare",action="store_true")
    ap.add_argument("--inputs",type=Path)
    ap.add_argument("--cache",type=Path)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--fit",action="store_true")
    args = ap.parse_args()
    if args.prepare:
        prepare_refs(args.reference)
    if args.inputs:
        compare(args.inputs,args.reference,args.cache,args.out,args.fit)

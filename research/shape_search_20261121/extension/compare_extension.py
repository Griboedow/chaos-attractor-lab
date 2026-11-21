"""Compare each archived follow-up geometry with 636 previous samples and 177 references."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from compare import descriptors,downsample,write
from index_campaign import RESULTS,read

PACKET=RESULTS/"extension"


def descriptor_file(path):
    with np.load(path,allow_pickle=False) as f:
        arr=f["points"]
        assert arr.ndim==2 and arr.shape[1]==3 and np.isfinite(arr).all(),path
        return len(arr),{k:v.tolist() for k,v in descriptors(downsample(arr)).items()}


def baseline(cache,out):
    rows=[]
    for key,rel in read(RESULTS/"comparison-inputs.json").items():
        _,fps=descriptor_file(cache/rel)
        rows.append(dict(name=key,descriptors=fps))
    assert len(rows)==636
    write(out,rows)


def records(cache):
    rows={}
    def walk(val):
        if isinstance(val,list):
            for child in val: walk(child)
        if isinstance(val,dict):
            path=val.get("cloud_relative")
            sha=val.get("cloud_sha256") or val.get("sha256") or val.get("file_sha256")
            if path is None and sha and isinstance(val.get("cloud"),str):
                raw=Path(val["cloud"])
                if "extension" in raw.parts:
                    path=Path(*raw.parts[raw.parts.index("extension"):]).as_posix()
            if path and sha:
                if path in rows: assert rows[path]==sha,path
                rows[path]=sha
            else:
                for child in val.values(): walk(child)
    for directory in ("brunnian","projection","flows","maps"):
        for entry in sorted((PACKET/directory).rglob("*.json")):
            walk(read(entry))
    assert all(str(k).startswith("extension/") and (cache/k).is_file() for k in rows)
    return rows


def nearest(book,mode,fp):
    names=[r["name"] for r in book if mode in r["descriptors"]]
    arr=np.array([r["descriptors"][mode] for r in book if mode in r["descriptors"]])
    delta=arr-np.asarray(fp)
    scores=np.linalg.norm(delta[:,:65],axis=1)+.5*np.linalg.norm(delta[:,65:68],axis=1)+.25*np.linalg.norm(delta[:,68:],axis=1)
    return [(names[i],float(scores[i])) for i in np.argsort(scores)[:3]]


def compare(cache,archive):
    old=read(archive)
    references=read(cache/"extended-reference.json")["rows"]
    assert len(references)==177
    receipt=[]
    for rel,sha in sorted(records(cache).items()):
        path=cache/rel
        with path.open("rb") as f: assert hashlib.file_digest(f,"sha256").hexdigest()==sha,rel
        n,modes=descriptor_file(path)
        results={mode:{"prior_campaign":nearest(old,mode,fp),"reference":nearest(references,mode,fp)}
                 for mode,fp in modes.items()}
        receipt.append(dict(cloud=rel,sha256=sha,points=n,compared=results))
    write(PACKET/"geometry-comparison.json.gz",receipt)
    totals={}
    for row in receipt:
        group=row["cloud"].split("/")[1]
        t=totals.setdefault(group,{"clouds":0,"points":0})
        t["clouds"]+=1;t["points"]+=row["points"]
    summary=dict(previous_samples=len(old),reference_samples=len(references),by_family=totals,
                 clouds=len(receipt),points=sum(r["points"] for r in receipt))
    write(PACKET/"geometry-totals.json",summary)
    print(summary)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",required=True,type=Path);ap.add_argument("--baseline",required=True,type=Path)
    ap.add_argument("--prepare",action="store_true")
    args=ap.parse_args()
    if args.prepare: baseline(args.cache,args.baseline)
    else: compare(args.cache,args.baseline)

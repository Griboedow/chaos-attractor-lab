"""Audit the point count, finiteness and independently stored hash of each saved geometry."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from index_campaign import read,write,RESULTS


def expected():
    result={}
    def add(rows,hash_key="sha256"):
        for row in rows:
            rel=row["cloud_relative"]
            sha=row[hash_key]
            if rel in result:
                assert result[rel]==sha,rel
            result[rel]=sha
    packet=read(RESULTS/"flows/manifest.json")
    for name in ("cases","stable_cycles","linked_cycles","knotted_periodic_skeletons"):
        add(packet[name])
    add(read(RESULTS/"flows/nested-skeletons.json")["cores"])
    add([read(RESULTS/"flows/final-link-crosscheck.json")["symmetric_companion"]])
    with (RESULTS/"flows/periodic-catalog.csv").open(encoding="utf-8",newline="") as f:
        add(csv.DictReader(f))
    add(read(RESULTS/"maps/manifest.json")["finalists"])
    for row in read(RESULTS/"projection/manifest.json"):
        path=RESULTS/"projection"/row["record"]
        data=read(path)
        result[row["cloud_relative"]]=data["cloud_sha256"]
    for family,name,data,key in [("braids","manifest.json","rows","file_sha256"),
                                 ("delay","manifest.json",None,"file_sha256"),
                                 ("delay","sink-clouds.json",None,"file_sha256"),
                                 ("pinched","manifest.json","selected","file_sha256")]:
        rows=read(RESULTS/family/name)
        add(rows if data is None else rows[data],key)
    return result


def audit(cache):
    inputs=read(RESULTS/"comparison-inputs.json")
    known=expected()
    for path in (RESULTS/"projection/three/companions").glob("companion-*.json"):
        row=read(path)
        known[Path(row["cloud"]).relative_to(cache).as_posix()]=row["cloud_sha256"]
    rows=[]
    for case,path in sorted(inputs.items()):
        raw=(cache/path).read_bytes()
        digest=hashlib.sha256(raw).hexdigest()
        assert known[path]==digest,(case,path)
        with np.load(io.BytesIO(raw),allow_pickle=False) as saved:
            z=saved["points"]
        assert z.ndim==2 and z.shape[1]==3 and len(z)>3 and np.isfinite(z).all(),case
        rows.append({"id":case,"file_sha256":digest,"points":len(z),"storage_bytes":len(raw),"dtype":str(z.dtype)})
    groups={}
    for row in rows:
        group=row["id"].split("/")[0]
        totals=groups.setdefault(group,{"clouds":0,"points":0,"storage_bytes":0})
        totals["clouds"]+=1
        totals["points"]+=row["points"]
        totals["storage_bytes"]+=row["storage_bytes"]
    write(RESULTS/"cloud-audit.json.gz",rows)
    write(RESULTS/"cloud-totals.json",{"clouds":len(rows),"points":sum(v["points"] for v in rows),
                       "storage_bytes":sum(v["storage_bytes"] for v in rows),"by_group":groups})
    print(len(rows),"finite three-coordinate clouds with matching source checksums",flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",required=True,type=Path)
    audit(ap.parse_args().cache)

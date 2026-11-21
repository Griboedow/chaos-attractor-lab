"""Verify local cloud integrity, evidence paths and archive coverage."""
import hashlib
import json
import re
from pathlib import Path
from search import CACHE,RESULTS
from validate import write_json
import numpy as np


def audit():
    status=json.loads((RESULTS/"manifest.json").read_text())
    ids=set()
    for item in status["finalists"]:
        case=item["case"]
        if case in ids:
            raise ValueError("Repeated case")
        ids.add(case)
        path=Path(item["cloud"])
        if path.resolve()!=(CACHE.parent/item["cloud_relative"]).resolve():
            raise ValueError(f"{case}: paths differ")
        if hashlib.sha256(path.read_bytes()).hexdigest()!=item["sha256"]:
            raise ValueError(f"{case}: cloud checksum mismatch")
        with np.load(path) as archive:
            pts=archive[item["key"]]
            if pts.shape!=(240000,3) or pts.dtype!=np.float32 or not np.isfinite(pts).all():
                raise ValueError(f"{case}: invalid point cloud")
        for name,value in item.items():
            if isinstance(value,str) and value.endswith((".json",".md")) and not (RESULTS/value).exists():
                raise FileNotFoundError(f"{case}: missing evidence {name} {value}")
    for name in ("geometry.json","all-archive-comparison.json"):
        other=json.loads((RESULTS/name).read_text())
        if {row["case"] for row in other["comparisons"]}!=ids:
            raise ValueError(f"{name}: different case set")
    for path in (RESULTS/"README.md",RESULTS/"dimension-proof.md"):
        for target in re.findall(r"\]\(([^)]+)\)",path.read_text(encoding="utf-8")):
            if "://" not in target and not (RESULTS/target).exists():
                raise FileNotFoundError(f"{path.name}: missing {target}")
    for path in Path(__file__).parent.glob("*.py"):
        compile(path.read_text(encoding="utf-8"),str(path),"exec")
    oversized=[p.name for p in RESULTS.rglob("*.json") if p.stat().st_size>=256*1024]
    if oversized:
        raise ValueError(f"Oversized result files: {oversized}")
    write_json(RESULTS/"packet-audit.json",dict(checked_clouds=len(ids),
        checked_points=len(ids)*240000,points_per_cloud=240000,coordinates=3,dtype="float32",
        all_finite=True,absolute_and_relative_paths_agree=True,all_sha256_match=True,
        all_linked_evidence_found=True,both_archive_comparisons_cover_all_cases=True,
        compiled_source_files=len(list(Path(__file__).parent.glob("*.py"))),
        maximum_tracked_json_bytes=max(p.stat().st_size for p in RESULTS.rglob("*.json"))))
    print("Audited",len(ids),"clouds",flush=True)


if __name__=="__main__":
    audit()

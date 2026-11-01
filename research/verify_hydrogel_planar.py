"""Verify the compact hydrogel packet; optionally repeat its linear controls."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKET = ROOT / "results/hydrogel-planar-controls"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--numerical",action="store_true")
    args = parser.parse_args()
    manifest = json.loads((PACKET/"manifest.json").read_text(encoding="utf-8"))
    for item in manifest["files"]:
        path = PACKET/item["path"]
        assert path.stat().st_size == item["bytes"], path
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"], path
    for item in manifest["dependencies"]:
        path = ROOT/item["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"], path
    runs = {p.stem:json.loads(p.read_text(encoding="utf-8")) for p in (PACKET/"runs").glob("*.json")}
    assert len(runs) == 15 and all(r["success"] for r in runs.values())
    arrays = json.loads((PACKET/"arrays.json").read_text(encoding="utf-8"))["files"]
    assert len(arrays) == len(runs)
    for item in arrays:
        r = runs[item["run"].replace("/","-")]
        assert item["arrays"]["y"]["shape"] == [2*r["segments"],round(50*r["duration"])+1]
    inspections = json.loads((PACKET/"inspection.json").read_text(encoding="utf-8"))
    scouts = [r for r in inspections if "scout-01" in r["folder"]]
    assert len(inspections) == 15 and len(scouts) == 8
    assert sum(r["samples_with_proper_crossing"] > 0 for r in scouts) == 6
    assert all(r["late_samples"] == 2001 for r in scouts)
    assert abs(runs["event-32-BDF"]["period_mean"]-runs["event-32-Radau"]["period_mean"]) < 2e-9
    for name,bound in [("negative-32-BDF",5.35e-7),("negative-64-Radau",2e-8),("long-127-64",1.62e-6)]:
        assert runs[name]["full_state_return_max_last10"] < bound
    numerical = "not requested"
    if args.numerical:
        import numpy as np
        from hydrogel_planar import linear_controls
        actual = linear_controls()
        saved = json.loads((PACKET/"linear.json").read_text(encoding="utf-8"))
        for a,b in zip(actual["rows"],saved["rows"],strict=True):
            assert (a["segments"],a["plane"]) == (b["segments"],b["plane"])
            for method in ("hermite","angle"):
                np.testing.assert_allclose(a["critical_positive_chi"][method],b["critical_positive_chi"][method],rtol=1e-6,atol=1e-6)
                np.testing.assert_allclose(a["chi100_"+method],b["chi100_"+method],rtol=1e-6,atol=1e-6)
        np.testing.assert_allclose(actual["passive"]["hermite"],saved["passive"]["hermite"],rtol=1e-6,atol=1e-8)
        numerical = "both spatial methods, passive control, RHS linearization and independent drag quadrature passed"
    print(json.dumps(dict(retained_hashes=len(manifest["files"])+len(manifest["dependencies"]),
                          compact_runs=len(runs),scouts=len(scouts),external_array_bytes=sum(a["bytes"] for a in arrays),
                          numerical=numerical,scope="No nonlinear replay, external-array verification or continuum certificate in this command.")))


if __name__ == "__main__":
    main()

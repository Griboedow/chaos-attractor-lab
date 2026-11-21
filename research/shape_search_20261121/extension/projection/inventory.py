"""Verify the external checkpoints and native clouds; publish portable compact indexes."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_text(encoding="utf8"))


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf8")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    output, cache = args.output, args.cache
    anchors = read(args.study / "three/manifest.json")
    labels = {409700: "Twin density plumes", 1403259: "Paired folded petals",
              2298665: "Lobes with intermittent phase slips", 2502898: "Sheared twin folds"}
    cases = sorted((read(p) for p in (output / "cases").glob("case-*.json")), key=lambda row: row["case"])
    stats = {}
    checks = []
    for title, screen, deep in (("affine", "screen", "deep"), ("sheared", "screen-warp", "deep-warp")):
        stats[title] = read(cache / deep / "stage.json")
        with np.load(cache / deep / "shortlist.npz") as shortlist:
            assert len(set(shortlist["ids"].tolist())) == stats[title]["shortlist"]
        rows = [row for row in cases if (row["case"] >= 2000000) == (title == "sheared")]
        assert len(rows) == stats[title]["shortlist"]
        stats[title]["cpu_records"] = len(rows)
        stats[title]["retained"] = sum(row["retained"] for row in rows)
        for path in sorted((cache / screen).glob("screen-*.npz")) + [cache / deep / "shortlist.npz", cache / deep / "stage.json"]:
            checks.append({"path": path.relative_to(cache).as_posix(), "sha256": digest(path)})
    kept = []
    for doc in cases:
        if not doc["retained"]:
            continue
        case = doc["case"]
        path = cache / "clouds" / f"affine-{case}.npz"
        assert digest(path) == doc["cloud_sha256"]
        with np.load(path) as native:
            assert native["points"].shape == (250000, 3) and native["points"].dtype == np.float32
            assert int(native["case"]) == case
        record = f"cases/case-{case}.json"
        more = {}
        for folder in ("precision", "precision-long", "extended"):
            found = list((output / folder).glob(f"*-{case}.json"))
            if found:
                more[folder] = found[0].relative_to(output).as_posix()
        kept.append(dict(case=case, family=doc["family"], label=labels.get(case, doc["family"] + f" {case}"),
                         study="extension/projection", cloud=path.resolve().as_posix(), cloud_relative=doc["cloud_relative"],
                         record=record, rates=doc["rates"], minimum=doc["minimum"], three_positive_starts=doc["cohort_three"],
                         total_starts=len(doc["cohort"]), observed=doc["windows"] * doc["window_length"],
                         points_saved=doc["points_saved"], burn=doc["burn"], cloud_sha256=doc["cloud_sha256"],
                         nearest_old_case=anchors[doc["nearest_old_index"]]["case"],
                         old_distance=doc["normalized_descriptor_distance_to_nearest_old_mean"], controls=more))
    write(output / "manifest.json", kept)
    write(output / "checkpoints.json", checks)
    write(output / "summary.json", dict(streams=stats, unique_new=sum(v["unique"] for v in stats.values()),
          retained=len(kept), cpu_rejected=[v["case"] for v in cases if not v["retained"]],
          checkpoint_index_sha256=digest(output / "checkpoints.json"), manifest_sha256=digest(output / "manifest.json")))
    print(stats)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for name in ("study", "cache", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    main(parser.parse_args())

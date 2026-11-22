"""Reconstruct the search counts, selected native clouds, and temporal diagnostics."""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from survey import features

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    m = np.load(a.cache / "screen.npz")
    ev, _ = features(m["params"], m["metrics"])
    positive = m["metrics"][:, :4].min(1) > 0.012
    deep = np.load(a.cache / "full-spectra.npz")
    out = dict(
        screened=len(positive),
        short_positive=int(positive.sum()),
        short_positive_and_covariance=int(np.sum(positive & (ev[:, 0] > 0.0003))),
        deep=len(deep["ids"]),
        deep_minimum_positive=int(len(np.load(a.cache / "morphology.npz")["ids"])),
        deep_two_rates_over_001=int(np.sum(deep["rates"][:, :, 1].min(1) > 0.01)),
        visual_shortlist=144,
        visual_rechecks=[
            24523,
            27148,
            86780,
            16926,
            1935,
            68935,
            45941,
            56729,
            129100,
            41561,
            75157,
            43025,
            97907,
            101804,
            30624,
            66075,
            58164,
            70380,
        ],
        records={},
    )
    for label in ["45941", "86780", "86780-upper", "56729"]:
        path = a.cache / f"case-{label}.npz"
        q = np.load(path)
        z = q["points"].astype("f8")
        row = dict(
            npz=path.name,
            points=len(z),
            parameters=q["params"].tolist(),
            initial=q["initial"].tolist(),
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        )
        if label.startswith("86780"):
            row["field_chord_medians"] = {
                str(k): float(np.median(np.linalg.norm(z[k:, :2] - z[:-k, :2], axis=1)))
                for k in [1, 2, 3, 6]
            }
            row["field_chord_fraction_under_002"] = {
                str(k): float(
                    np.mean(np.linalg.norm(z[k:, :2] - z[:-k, :2], axis=1) < 0.02)
                )
                for k in [1, 2, 3, 6]
            }
        if label == "56729":
            row["real_field_ranges_by_parity"] = [
                [float(z[i::2, 0].min()), float(z[i::2, 0].max())] for i in range(2)
            ]
        out["records"][label] = row
    a.out.write_text(json.dumps(out, indent=2) + "\n")

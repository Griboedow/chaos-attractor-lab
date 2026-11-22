"""Freeze or reproduce a trajectory; compare long windows and independent starts."""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from scipy.spatial.distance import jensenshannon
from model import trace, independent

NAMES = (
    "pump",
    "field_retention",
    "base_detuning",
    "fast_shift",
    "memory_gain",
    "relaxation",
)
INITIAL = [0.173, -0.341, 0.271]


def distribution(z, reference):
    edges = [
        np.linspace(a, b, 65)
        for a, b in zip(reference.min(0) - 1e-5, reference.max(0) + 1e-5)
    ]

    def histogram(x):
        return [
            np.histogram2d(
                np.clip(x[:, i], edges[i][0], edges[i][-1]),
                np.clip(x[:, j], edges[j][0], edges[j][-1]),
                bins=[edges[i], edges[j]],
            )[0].ravel()
            for i, j in [(0, 1), (0, 2), (1, 2)]
        ]

    return [
        float(jensenshannon(a, b, base=2))
        for a, b in zip(histogram(z), histogram(reference))
    ]


def diagnostics(z, case):
    if case == 45941:
        specs = [
            dict(height=h, distance=10, prominence=0.1) for h in [-1.55, -1.5, -1.45]
        ]
    elif case == 86780:
        specs = [
            dict(height=h, distance=180, prominence=0.25) for h in [5.64, 5.66, 5.68]
        ]
    else:
        return {}
    ret = []
    for spec in specs:
        peaks = find_peaks(z[:, 2], **spec)[0]
        lags = np.diff(peaks)
        ret.append(
            dict(
                rule=spec,
                events=len(peaks),
                interval_quantiles=np.quantile(
                    lags, [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
                ).tolist()
                if len(lags)
                else [],
                interval_counts={
                    str(int(i)): int(n)
                    for i, n in zip(*np.unique(lags, return_counts=True))
                    if n >= 0.01 * len(lags)
                },
            )
        )
    return dict(thermal_peak_intervals=ret)


def run(a):
    a.cache.mkdir(parents=True, exist_ok=True)
    if a.record:
        rec = json.loads(a.record.read_text())
        case = rec["case"]
        params = np.array([rec["parameters"][x] for x in NAMES])
        init = np.array(rec["initial"])
    else:
        case = a.case
        params = np.load(a.cache / "screen.npz")["params"][case].astype("f8")
        init = np.array(INITIAL)
    burn = 131072
    window = 262144
    windows = 12
    keep = 786432
    z, rates = trace(params, init, burn, window, windows, keep)
    cloud = a.cache / f"case-{case}.npz"
    np.savez_compressed(cloud, points=z.astype("f4"), params=params, initial=init)
    if a.cloud_only:
        print(cloud, hashlib.sha256(cloud.read_bytes()).hexdigest())
        return
    rec = dict(
        case=case,
        parameters=dict(zip(NAMES, params.tolist())),
        initial=init.tolist(),
        burn=burn,
        window=window,
        windows=windows,
        kept=keep,
        validated_orbit_steps=window * windows,
        xyz=["Re field", "Im field", "dynamic detuning"],
        npz=cloud.name,
        sha256=hashlib.sha256(cloud.read_bytes()).hexdigest(),
        rates=rates.tolist(),
        mean_rates=rates.mean(0).tolist(),
        analytic_log_determinant=float(np.log(params[1] ** 2 * (1 - params[5]))),
        max_volume_discrepancy=float(
            np.max(abs(rates.sum(1) - np.log(params[1] ** 2 * (1 - params[5]))))
        ),
        bounds=np.stack([z.min(0), z.max(0)]).tolist(),
        spatial_mean=z.mean(0).tolist(),
        spatial_std=z.std(0).tolist(),
        diagnostics=diagnostics(z, case),
    )
    rng = np.random.default_rng(229866512)
    starts = [init] + [
        np.r_[
            rng.uniform(-params[0] / (1 - params[1]), params[0] / (1 - params[1]), 2),
            rng.uniform(min(0, params[4]), max(0, params[4])),
        ]
        for _ in range(16)
    ]
    rec["start_sampling"] = (
        "PCG64(229866512); x,y independent uniform in [-P/(1-a),P/(1-a)], T uniform between 0 and g."
    )
    rec["other_start_controls"] = []
    for start in starts:
        orb, rr = trace(params, start, 131072, 32768, 4, 65536)
        ref = z[-65536:]
        rec["other_start_controls"].append(
            dict(
                initial=start.tolist(),
                spectrum=rr.mean(0).tolist(),
                distances=distribution(orb, ref),
                means=orb.mean(0).tolist(),
            )
        )
    y = independent(params, init, 131072, 786432, 786432)
    rec["complex_arithmetic"] = dict(
        n=len(y),
        distances=distribution(y, z),
        mean=y.mean(0).tolist(),
        std=y.std(0).tolist(),
        diagnostics=diagnostics(y, case),
    )
    rec["split_reference_distances"] = distribution(z[: len(z) // 2], z[len(z) // 2 :])
    frozen = float(z[:, 2].mean())
    _, cf = trace(params, init, 131072, 65536, 4, 0, frozen)
    rec["fixed_detuning_counterfactual"] = dict(
        T=frozen, mean_rates=cf.mean(0).tolist()
    )
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / f"case-{case}.json").write_text(json.dumps(rec, indent=2) + "\n")
    print(
        case,
        rec["mean_rates"],
        rec["complex_arithmetic"]["distances"],
        [k["spectrum"][0] for k in rec["other_start_controls"]],
        flush=True,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("."))
    ap.add_argument("--case", type=int)
    ap.add_argument("--record", type=Path)
    ap.add_argument("--cloud-only", action="store_true")
    run(ap.parse_args())

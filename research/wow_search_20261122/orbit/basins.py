"""Test separation, sample the initial-condition plane, and reproduce the upper branch."""

import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
from numba import njit, prange, set_num_threads
from model import trace, independent
from validate import distribution, NAMES


@njit(cache=True, parallel=True)
def scan(params, starts, burn=32768, steps=8192):
    out = np.empty((len(starts), 4))
    for i in prange(len(starts)):
        z = independent(params, starts[i], burn, steps, steps)
        out[i] = np.array([z[:, 2].min(), z[:, 2].mean(), z[:, 2].max(), z[:, 2].std()])
    return out


@njit(cache=True, parallel=True)
def long_blocks(params, starts, n=16384, blocks=1024):
    out = np.empty((len(starts), blocks, 4))
    for k in prange(len(starts)):
        x = independent(params, starts[k], 131072, 1, 1)[-1]
        for b in range(blocks):
            z = independent(params, x, 0, n, n)
            x = z[-1]
            out[k, b] = np.array(
                [z[:, 2].min(), z[:, 2].mean(), z[:, 2].max(), z[:, 2].std()]
            )
    return out


def main(a):
    a.cache.mkdir(parents=True, exist_ok=True)
    rec = json.loads(a.record.read_text())
    params = np.array([rec["parameters"][k] for k in NAMES])
    set_num_threads(a.threads)
    low = np.array(rec["initial"])
    high = np.array(rec["other_start_controls"][1]["initial"])
    n = 131072
    w = 262144
    nr = 12
    keep = 786432
    z, rr = trace(params, high, n, w, nr, keep)
    target = a.cache / "case-86780-upper.npz"
    np.savez_compressed(target, points=z.astype("f4"), params=params, initial=high)
    if a.cloud_only:
        print(target, hashlib.sha256(target.read_bytes()).hexdigest())
        return
    ind = independent(params, high, n, keep, keep)
    out = dict(
        initial_upper=high.tolist(),
        lower_npz="case-86780.npz",
        upper_npz=target.name,
        upper_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        upper_rates=rr.tolist(),
        upper_bounds=np.stack([z.min(0), z.max(0)]).tolist(),
        independent_upper_distances=distribution(ind, z),
        independent_upper_mean=ind.mean(0).tolist(),
        lower_bounds=rec["bounds"],
    )
    gen = np.random.default_rng(229866513)
    xy = np.vstack(
        [
            low[:2],
            np.zeros(2),
            gen.uniform(
                -params[0] / (1 - params[1]), params[0] / (1 - params[1]), (7, 2)
            ),
        ]
    )
    thermal = np.linspace(0, params[4], 257)
    starts = np.array([[x, y, t] for x, y in xy for t in thermal])
    ss = scan(params, starts).reshape((len(xy), len(thermal), 4))
    threshold = 7.5
    lab = np.where(ss[:, :, 2] < threshold, 0, np.where(ss[:, :, 0] > threshold, 1, 2))
    out["basin_scan"] = dict(
        field_initials=xy.tolist(),
        thermal_min=0,
        thermal_max=params[4],
        thermal_grid_size=len(thermal),
        burn=32768,
        steps=8192,
        separating_threshold=threshold,
        labels="0 maximum T < 7.5; 1 minimum T > 7.5; 2 window straddles 7.5",
        observed_counts=np.bincount(lab.ravel(), minlength=3).tolist(),
        rows=["".join(map(str, row.tolist())) for row in lab],
        maximum_low=float(np.max(ss[:, :, 2][lab == 0])),
        minimum_high=float(np.min(ss[:, :, 0][lab == 1])),
    )
    seeds = np.vstack(
        [
            low,
            high,
            np.array(
                [rec["other_start_controls"][k]["initial"] for k in [2, 3, 5, 6, 9, 15]]
            ),
        ]
    )
    blocks = long_blocks(params, seeds)
    out["long_blocks"] = dict(
        test="independent complex implementation",
        block_size=16384,
        blocks_per_start=1024,
        transient=131072,
        starts=[
            dict(
                initial=st.tolist(),
                min=float(v[:, 0].min()),
                max=float(v[:, 2].max()),
                range_of_block_mean=[float(v[:, 1].min()), float(v[:, 1].max())],
                crossing_blocks=int(
                    np.sum((v[:, 0] < threshold) & (v[:, 2] > threshold))
                ),
            )
            for st, v in zip(seeds, blocks)
        ],
    )
    np.savez_compressed(
        a.cache / "case-86780-basins.npz",
        starts=starts,
        statistics=ss,
        block_starts=seeds,
        blocks=blocks,
    )
    (a.out / "case-86780-basins.json").write_text(json.dumps(out, indent=2) + "\n")
    print(out["basin_scan"]["observed_counts"], out["long_blocks"], flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--record", required=True, type=Path)
    p.add_argument("--cache", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--threads", type=int, default=5)
    p.add_argument("--cloud-only", action="store_true")
    main(p.parse_args())

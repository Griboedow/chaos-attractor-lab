"""Target the order-three and order-five windows of the scalar cubic radial map."""

import numpy as np
from scipy.stats import qmc
from numba import set_num_threads
from icons import CACHE
from nonpower import candidate

BASE = 200000
POWER = 13
SEED = 2026112231


def parameters_odd():
    u = qmc.Sobol(8, scramble=True, seed=SEED).random_base2(POWER)
    rows = np.empty((len(u), 9))
    centers = [2.4531555836, 2.3943931930, 2.5059635817]
    for i, t in enumerate(u):
        n = 8 + int(15 * t[0])
        l = centers[int(3 * t[1])] + 0.022 * (t[2] - 0.27)
        a = -l * (1.12 + 0.24 * t[3]) ** 2
        rows[i] = (
            n,
            l,
            a,
            0.24 * (t[4] - 0.5),
            0.04 + 0.92 * t[5],
            0.22 * (t[6] - 0.5),
            0.25 * (t[7] - 0.5),
            0.27 * np.cos(2 * np.pi * t[4]),
            0.27 * np.sin(2 * np.pi * t[4]),
        )
    return rows


if __name__ == "__main__":
    set_num_threads(3)
    rows = parameters_odd()
    chunk = 1024
    metrics = []
    for a in range(0, len(rows), chunk):
        ix = np.arange(a, a + chunk)
        metrics.append(candidate(rows, ix))
    s = np.vstack(metrics)
    good = np.flatnonzero(s[:, 0])
    print("passed", len(good), flush=True)
    rank = sorted(good, key=lambda k: -s[k, 2] * min(1.0, s[k, 1] / 0.018))
    selected = np.r_[rank[:25], [k for k in rank if s[k, 0] == 5][:8]].astype(int)
    print([(int(k + BASE), s[k].round(5).tolist()) for k in selected[:16]], flush=True)
    np.savez_compressed(
        CACHE / "odd-scout.npz", parameters=rows, metrics=s, selected=selected
    )

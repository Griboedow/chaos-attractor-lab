"""Find separate radius intervals for cyclic phases of retained scout rows."""

import numpy as np
from icons import CACHE, INITIAL, orbit
from numba import njit, prange, set_num_threads


@njit(parallel=True)
def candidate(rows, ids):
    result = np.zeros((len(ids), 5))
    periods = (3, 5, 6, 7, 9, 10, 12)
    for j in prange(len(ids)):
        p = rows[ids[j]]
        z, rate, ok = orbit(p, INITIAL, 10000, 12600, 6)
        if not ok or rate[:, 0].min() < 0.001:
            continue
        rad = np.abs(z)
        best = 0.0
        for q in periods:
            low = np.array([rad[k::q].min() for k in range(q)])
            high = np.array([rad[k::q].max() for k in range(q)])
            idx = np.argsort(low)
            gap = np.min(low[idx[1:]] - high[idx[:-1]])
            if gap > best:
                best = gap
                result[j] = np.array([q, gap, rate[:, 0].mean(), rad.min(), rad.max()])
    return result


if __name__ == "__main__":
    set_num_threads(3)
    data = np.load(CACHE / "scout-14.npz")
    ids = np.flatnonzero(data["metrics"][:, 0])
    s = candidate(data["parameters"], ids)
    print("nonpowers strict sampled bands", np.count_nonzero(s[:, 0]), flush=True)
    print(
        [
            (int(ids[k]), s[k].round(6).tolist())
            for k in np.argsort(s[:, 1])[::-1][:18]
            if s[k, 1] > 0
        ],
        flush=True,
    )
    np.savez_compressed(CACHE / "nonpower.npz", ids=ids, metrics=s)

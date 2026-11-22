"""Compute two-fold-orthogonalized tangent spectra across the deep shortlist."""

import argparse
from pathlib import Path
import numpy as np
from numba import njit, prange, set_num_threads
from model import trace


@njit(parallel=True, cache=True)
def spectra(params, initial):
    out = np.empty((len(params), 4, 3))
    for i in prange(len(params)):
        _, rates = trace(params[i], initial, 8192, 4096, 4, 0)
        out[i] = rates
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--threads", type=int, default=4)
    a = p.parse_args()
    set_num_threads(a.threads)
    q = np.load(a.cache / "deep-design.npz")
    params = q["params"].astype("f8")
    scores = spectra(params, np.array([0.173, -0.341, 0.271]))
    np.savez_compressed(
        a.cache / "full-spectra.npz", ids=q["ids"], params=params, rates=scores
    )
    i = np.argsort(scores[:, :, 1].min(1))[::-1][:12]
    print(list(zip(q["ids"][i].tolist(), scores[i].mean(1).tolist())), flush=True)

"""Store the exact consecutive binary64 states of Case 2298665 externally."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
import sys

import numpy as np
from numba import njit

EXTENSION = Path(__file__).resolve().parents[2] / "extension" / "projection"
sys.path.insert(0, str(EXTENSION))
from affine_model import step
from confirm import from_record


@njit(cache=True)
def advance(raw, state, length, output, save):
    for k in range(length):
        state, _ = step(state, raw)
        if save:
            output[k, :] = state
    return state


def initial_state(doc, cohort):
    if cohort < len(doc["cohort"]):
        return np.array(doc["cohort"][cohort]["initial"])
    seed = np.random.SeedSequence([472919, 2298665, cohort])
    return np.random.default_rng(seed).uniform(-2, 2, 3)


def simulate(args):
    record, folder, cohort, length, burn = args
    path = folder / f"trace-{cohort:02d}-{length}.npy"
    if path.exists():
        raise FileExistsError(path)
    work = path.with_suffix(".partial.npy")
    doc = json.loads(record.read_text())
    raw = from_record(doc)
    init = initial_state(doc, cohort)
    state = advance(raw, init, burn, np.empty((0, 3)), False)
    dest = np.lib.format.open_memmap(work, mode="w+", dtype="float64", shape=(length, 3))
    for start in range(0, length, 1000000):
        block = dest[start:start + 1000000]
        state = advance(raw, state, len(block), block, True)
    dest.flush()
    del block, dest
    os.replace(work, path)
    print(path, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--record", required=True, type=Path)
    p.add_argument("--cache", required=True, type=Path)
    p.add_argument("--cohorts", nargs="+", type=int, required=True)
    p.add_argument("--length", type=int, required=True)
    p.add_argument("--burn", type=int, default=100000)
    p.add_argument("--workers", type=int, default=1)
    a = p.parse_args()
    a.cache.mkdir(parents=True, exist_ok=True)
    calls = [(a.record, a.cache, i, a.length, a.burn) for i in a.cohorts]
    if a.workers == 1:
        for call in calls:
            simulate(call)
    else:
        with ProcessPoolExecutor(a.workers) as pool:
            list(pool.map(simulate, calls))

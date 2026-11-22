"""Classify successive section visits relative to eight orientations of the motor braid."""

import argparse
import json
from itertools import permutations
from pathlib import Path

import numpy as np
from motor_rows import source_force_constraint_field
from numba import njit
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
TRAIN = ROOT / "results/motor-four-geometry/long/protocol.json"
LARGE = ROOT / "results/motor-four-b-study/motor-four-symmetric-long/protocol.json"
LEVELS = (-0.125, 0.0, 0.125)


@njit
def field(x, p, out):
    rho, r, b = p
    mean_s = (x[4] + x[5] + x[6] + x[7]) / 4
    for j in range(4):
        displacement = x[8 + j] if j < 3 else -x[8] - x[9] - x[10]
        v = -rho * displacement + r * (x[4 + j] - mean_s)
        out[j] = 1 - x[j] - v * x[4 + j]
        out[4 + j] = b - x[4 + j] + v * x[j]
        if j < 3:
            out[8 + j] = v


@njit
def sample(initial, p, warm, duration, step):
    stride = round(0.1 / step)
    assert abs(stride * step - 0.1) < 1e-12
    n_warm = round(warm / step)
    n_observe = round(duration / step)
    x = initial.copy()
    out = np.empty((n_observe // stride + 1, len(x)))
    k = np.empty((4, len(x)))
    stage = np.empty_like(x)
    for n in range(n_warm + n_observe + 1):
        if n >= n_warm and (n - n_warm) % stride == 0:
            out[(n - n_warm) // stride] = x
        if n == n_warm + n_observe:
            break
        field(x, p, k[0])
        for i in range(1, 4):
            width = step if i == 3 else step * 0.5
            for j in range(len(x)):
                stage[j] = x[j] + width * k[i - 1, j]
            field(stage, p, k[i])
        for j in range(len(x)):
            x[j] += step / 6 * (k[0, j] + 2 * k[1, j] + 2 * k[2, j] + k[3, j])
    return out


def section(states, p, level, spacing):
    """Upward mean-c crossings with cubic Hermite interpolation in full state."""
    q = states[:, :4].mean(axis=1)
    ids = np.flatnonzero((q[:-1] < level) & (q[1:] >= level))
    x, y = states[ids], states[ids + 1]
    dx, dy = np.empty_like(x), np.empty_like(y)
    for j in range(len(x)):
        field(x[j], p, dx[j])
        field(y[j], p, dy[j])
    c1 = spacing * dx
    c2 = 3 * (y - x) - spacing * (2 * dx + dy)
    c3 = 2 * (x - y) + spacing * (dx + dy)
    coefficients = np.stack([c[:, :4].mean(axis=1) for c in (x, c1, c2, c3)])
    lo, hi = np.zeros(len(x)), np.ones(len(x))
    for _ in range(32):
        mid = (lo + hi) / 2
        above = (
            (coefficients[3] * mid + coefficients[2]) * mid + coefficients[1]
        ) * mid + coefficients[0] >= level
        hi = np.where(above, mid, hi)
        lo = np.where(above, lo, mid)
    t = ((lo + hi) / 2)[:, None]
    event = ((c3 * t + c2) * t + c1) * t + x
    assert np.max(np.abs(event[:, :4].mean(axis=1) - level)) < 1e-9
    full = np.column_stack([event, -event[:, 8:].sum(axis=1)])
    weights = np.r_[np.ones(8), np.full(4, p[0] / p[1])] / np.sqrt(1 + 1.5**2)
    return (ids + t[:, 0]) * spacing, full * weights


def templates(points, row):
    groups = [[] for _ in range(8)]
    for perm in permutations(range(4)):
        parity = sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4)) % 2
        key = 2 * perm.index(row) + parity
        indices = np.array(perm)
        groups[key].append(points[:, np.r_[indices, 4 + indices, 8 + indices]])
    return [cKDTree(np.concatenate(group)) for group in groups]


def distances(trees, points):
    return np.column_stack([tree.query(points)[0] for tree in trees])


def visits(d, time, radius, ratio):
    ranks = np.argsort(d, axis=1)
    best, second = np.take_along_axis(d, ranks[:, :2], axis=1).T
    take = (best < radius) & (second > ratio * best)
    labels = ranks[take, 0]
    before, after = labels[:-1], labels[1:]
    change = before != after
    a, b = before[change], after[change]
    edge = (a // 2 != b // 2) & (a % 2 != b % 2)
    face = (a // 2 != b // 2) & (a % 2 == b % 2)
    opposite = a // 2 == b // 2
    matrix = np.zeros((8, 8), dtype=int)
    np.add.at(matrix, (a, b), 1)
    return {
        "section_count": len(d),
        "visits": np.bincount(labels, minlength=8).tolist(),
        "changes": len(a),
        "edge": int(edge.sum()),
        "face": int(face.sum()),
        "opposite": int(opposite.sum()),
        "median_elapsed": float(np.median(np.diff(time[take])[change]))
        if len(a)
        else None,
        "matrix": matrix.tolist(),
    }


def simple(result):
    return {
        k: v
        for k, v in result.items()
        if k not in ("matrix", "median_elapsed", "visits")
    }


def settings(level):
    yield 0.075, 2.0
    if level == 0:
        for value in (1.5, 2.5, 3.0):
            yield 0.075, value
        yield 0.05, 2.0
        yield 0.1, 2.0


def run(output, step, archived=None):
    jobs = [
        next(job for job in json.loads(TRAIN.read_text())["jobs"] if job["id"] == 3142)
    ]
    jobs += json.loads(LARGE.read_text())["jobs"]
    tables = []
    for idx, job in enumerate(jobs):
        p = np.array(job["parameters"])
        x0 = np.array(job["initial"])
        derivative = np.empty(11)
        field(x0, p, derivative)
        np.testing.assert_allclose(
            derivative, source_force_constraint_field(x0, p, 4), rtol=1e-12, atol=1e-12
        )
        if archived is None:
            states = sample(x0, p, job["warm"], job["duration"], step)
            spacing = 0.1
        else:
            source = "motor-four-symmetric-long" if idx else "motor-four-long"
            with np.load(archived / source / f"case-{job['id']:04d}.npz") as data:
                np.testing.assert_allclose(p, data["parameters"], rtol=0, atol=0)
                states = data["states"]
                spacing = float(np.median(np.diff(data["time"][:1000])))
        if not idx:
            half = len(states) // 2
            x = states[:half, 8:]
            row = int(np.argmin(np.r_[x.mean(axis=0), -x.sum(axis=1).mean()]))
            trees = {
                level: templates(section(states[:half], p, level, spacing)[1], row)
                for level in LEVELS
            }
            states = states[half:]
        record = {"b": float(p[2]), "runs": []}
        for level in LEVELS:
            time, points = section(states, p, level, spacing)
            d = distances(trees[level], points)
            for radius, ratio in settings(level):
                v = visits(d, time, radius, ratio)
                primary = level == 0 and radius == 0.075 and ratio == 2
                entry = dict(
                    level=level,
                    radius=radius,
                    ratio=ratio,
                    **(v if primary else simple(v)),
                )
                if primary:
                    entry["halves"] = [
                        simple(visits(group_d, group_t, radius, ratio))
                        for group_d, group_t in zip(
                            np.array_split(d, 2), np.array_split(time, 2)
                        )
                    ]
                record["runs"].append(entry)
        tables.append(record)
        main = next(row for row in record["runs"] if "matrix" in row)
        print(
            p[2],
            "visits",
            sum(main["visits"]),
            "changes",
            main["changes"],
            "edges",
            main["edge"],
            flush=True,
        )
    result = {
        "method": "archived DOP853" if archived else "RK4",
        "step": None if archived else step,
        "original_distinguished_row": row,
        "training_duration": 6000,
        "wide_duration": 60000,
        "sampling": "Hermite-interpolated upward mean-c crossings",
        "class_order": "2*distinguished row + parity of permutation relative to the training history",
        "runs": tables,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=float, choices=(0.01, 0.005), default=0.01)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--archived",
        type=Path,
        help="Optional root containing original DOP853 full histories",
    )
    args = parser.parse_args()
    run(args.output, args.step, args.archived)

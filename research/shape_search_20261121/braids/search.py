"""Search Fourier braids in a contracting disk bundle over an expanding circle."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def terms(case):
    return [(k, complex(x, y)) for k, x, y in case["harmonics"]]


def center(case, angle):
    return sum(c * np.exp(1j * k * angle) for k, c in terms(case)) / case["scale"]


def draw_case(rng, case_id):
    degree = int(rng.choice([2, 3, 3, 4, 4, 5, 6, 7]))
    modes = np.r_[np.arange(-11, 0), np.arange(1, 12)]
    size = int(rng.integers(2, 6))
    powers = rng.choice(modes, size=size, replace=False)
    coeff = rng.integers(-8, 9, size=(size, 2))
    coeff[0] *= 2
    coeff[~np.any(coeff, axis=1), 0] = 1
    return {"case": case_id, "degree": degree,
            "harmonics": [[int(k), int(x), int(y)] for k, (x,y) in zip(powers, coeff)],
            "scale": 1, "base_ripple": 0.14, "base_offset": 0.137,
            "fiber_turns": int(rng.integers(-2, 3)), "fiber_phase": 0.31}


def clearance(case, nodes):
    """Grid-to-continuum estimates using a global bound on the Fourier derivative."""
    p = case["degree"]
    count = nodes * p
    grid = 2 * np.pi * np.arange(count) / count
    values = center(case, grid)
    speed_bound = sum(abs(k * c) for k, c in terms(case)) / case["scale"]
    interpolation = speed_bound * np.pi / count
    radius_upper = np.max(abs(values)) + interpolation
    gap_lower = min(np.min(abs(values - np.roll(values, shift * nodes)))
                    for shift in range(1, p // 2 + 1)) - 2 * interpolation
    return float(radius_upper), float(gap_lower)


def integer_determinant(matrix):
    """Fraction-free Gaussian elimination with row pivoting."""
    a = [list(map(int, row)) for row in matrix]
    previous, sign = 1, 1
    for k in range(len(a) - 1):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                numerator = a[k][k] * a[i][j] - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        previous = a[k][k]
    return sign * a[-1][-1]


def color_determinant(p, word):
    matrix = np.eye(p, dtype=object)
    for generator in word:
        i = abs(generator) - 1
        left, right = matrix[i].copy(), matrix[i + 1].copy()
        if generator > 0:
            matrix[i], matrix[i + 1] = 2 * left - right, left
        else:
            matrix[i], matrix[i + 1] = right, 2 * right - left
    relation = matrix - np.eye(p, dtype=object)
    return abs(integer_determinant(relation[:-1, :-1]))


def braid(case, nodes, camera=0.183, start=0.1731):
    p = case["degree"]
    angles = ((np.arange(nodes + 1)[:, None] / nodes + start) / p
              + np.arange(p)[None, :] / p) * (2 * np.pi)
    points = center(case, angles) * np.exp(1j * camera)
    order = points.real.argsort(axis=1)
    changes = np.flatnonzero(np.any(order[1:] != order[:-1], axis=1))
    word = []
    for k in changes:
        slots = np.flatnonzero(order[k] != order[k + 1])
        if len(slots) != 2 or slots[1] != slots[0] + 1:
            return None
        pos = slots[0]
        left, right = order[k, pos:pos + 2]
        dz = points[k:k + 2, left] - points[k:k + 2, right]
        mix = -dz[0].real / (dz[1].real - dz[0].real)
        height = dz[0].imag * (1 - mix) + dz[1].imag * mix
        word.append(int((pos + 1) * np.sign(height)))
    return word


def screen(count, seed, out):
    rng = np.random.default_rng(seed)
    survivors, keep, failures = [], {}, 0
    for index in range(count):
        case = draw_case(rng, index + 1)
        upper, lower = clearance(case, 96)
        case["scale"] = math.ceil(upper / 0.78)
        gap = lower / case["scale"]
        if gap < 0.155:
            continue
        case["fiber_radius"] = math.floor(1000 * min(0.19, 0.39 * gap)) / 1000
        word = braid(case, 512)
        if word is None:
            word = braid(case, 4096)
        if word is None:
            failures += 1
            continue
        determinant = color_determinant(case["degree"], word)
        assert determinant % 2 == 1
        case.update(core_determinant=determinant, crossings=len(word),
                    branch_gap_screen=gap, radius_upper_screen=upper / case["scale"])
        survivors.append(case)
        category = (case["degree"], determinant)
        if category not in keep or keep[category]["branch_gap_screen"] < gap:
            keep[category] = case
    write_json(out / "screen.json", {"draws": count, "seed": seed,
               "unresolved_crossings": failures, "survivors": len(survivors),
               "classes": len(keep), "leaders": list(keep.values())})
    write_json(out / "survivors.json", survivors)
    print(json.dumps({"draws": count, "survivors": len(survivors), "classes": len(keep)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=8000)
    parser.add_argument("--seed", type=int, default=20261121)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert color_determinant(2, [1] * 3) == 3
    assert color_determinant(3, [1, -2] * 2) == 5
    screen(args.count, args.seed, args.out)

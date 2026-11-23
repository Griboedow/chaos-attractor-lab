"""Interval shooting for the two period-four satellites of the driven cubic."""

import argparse
import itertools
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from elliptic_sites import OUT
from flint import acb, acb_poly, arb, arb_mat, arb_poly, ctx
from periodic_boxes import (
    cell,
    coefficient_norm,
    conjugate,
    drive_series,
    local_polynomials,
)
from scipy.integrate import solve_ivp

GAIN = "0.14260309809624272"
PANELS = 32
SUBSTEPS = 16
ORDER = 34
BITS = 208
RADIUS = 76


def proposal(row):
    k = float(GAIN)

    def f(t, h):
        b = 0.25 * np.sin(t) - 4 * np.sin(3 * t) + 3j * np.sin(6 * t)
        return -k * (6 * h + 3 * h**2).conj() * (b - 2 + 3 * h**2 + h**3)

    start = complex(*row["physical_offsets"][0])
    T = 2 * np.pi
    sol = solve_ivp(
        f,
        (0, 4 * T),
        [start],
        t_eval=np.arange(PANELS) * 4 * T / PANELS,
        method="DOP853",
        rtol=3e-14,
        atol=3e-25,
    )
    if not sol.success:
        raise ArithmeticError(sol.message)
    return [acb(1) + acb(float(z.real), float(z.imag)) for z in sol.y[0]]


def follow(nodes, cache, gain, step, radius, keep=False):
    n = len(nodes)
    substeps = 4 * len(cache) // n
    m = arb_mat(2 * n, 2 * n)
    residual = arb_mat(2 * n, 1)
    deriv = []
    boxes = []
    for i, z in enumerate(nodes):
        cols = [acb(1), acb(0, 1)]
        err = (arb(0), arb(2).sqrt() * radius, [arb(0), arb(0)])
        if keep:
            boxes.append([])
        for j in range(substeps):
            if keep:
                boxes[-1].append((z, err[1]))
            z, cols, err = cell(
                z, cols, err, cache[(i * substeps + j) % len(cache)], step, gain
            )
        target = nodes[(i + 1) % n]
        difference = target - z
        residual[2 * i, 0] = difference.real + arb(0, err[0])
        residual[2 * i + 1, 0] = difference.imag + arb(0, err[0])
        block = arb_mat(2, 2)
        for j, v in enumerate(cols):
            for k, value in enumerate((v.real, v.imag)):
                block[k, j] = value + arb(0, err[2][j])
                m[2 * i + k, 2 * i + j] = -block[k, j]
        m[2 * i, 2 * ((i + 1) % n)] += 1
        m[2 * i + 1, 2 * ((i + 1) % n) + 1] += 1
        deriv.append(block)
    return residual, m, deriv, boxes


def balanced_product(blocks):
    if len(blocks) == 1:
        return blocks[0]
    h = len(blocks) // 2
    return balanced_product(blocks[h:]) * balanced_product(blocks[:h])


def tube_bound(boxes, cache, k, h):
    """Enclose each whole arc in the same untwisted Floquet tube."""
    # Local variational columns are analytic in normalized time on |u|<=8.
    # Each real/imaginary entry has the displayed uniform Cauchy tail.
    size = 8 * h
    major = abs(k) * (
        12 + arb(3) / 2 * size.cosh() + 24 * (3 * size).cosh() + 18 * (6 * size).cosh()
    )
    tail = ((size * major).exp() * arb(8) ** (-ORDER) / 7).upper()

    def identity():
        return arb_mat([[1, 0], [0, 1]])

    def adj(m):
        return arb_mat([[m[1, 1], -m[0, 1]], [-m[1, 0], m[0, 0]]])

    core = [identity()]
    arcs = []
    for b, _ in cache:
        anti = [-6 * k * (v - (2 if j == 0 else 0)) for j, v in enumerate(b)]
        cols = []
        for start in (acb(1), acb(0, 1)):
            seq = [start]
            for n in range(ORDER):
                seq.append(
                    h
                    * sum(anti[j] * seq[n - j].conjugate() for j in range(n + 1))
                    / (n + 1)
                )
            cols.extend(
                [arb_poly([getattr(x, part) for x in seq]) for part in ("real", "imag")]
            )
        a, c, bb, dd = cols
        arcs.append((dd, -bb, -c, a))
        core.append(
            arb_mat(2, 2, [x(1) + arb(0, tail) for x in (a, bb, c, dd)]) * core[-1]
        )
    angle = arb(0)
    if not abs(k) * arb("55.5") * h < arb.pi() / 2:
        raise ArithmeticError("The core angle could wrap between samples")
    for start, end in itertools.pairwise(core):
        x, y = start[0, 0], start[1, 0]
        X, Y = end[0, 0], end[1, 0]
        dot = x * X + y * Y
        if not dot > 0:
            raise ArithmeticError("The core angle branch is unresolved")
        angle += ((x * Y - y * X) / dot).atan()
    core_turn = angle / (2 * arb.pi())
    if not arb(-3) / 2 < core_turn < arb(-1):
        raise ArithmeticError("The selected frame may twist around the phase tube")
    m = core[-1]
    co = (m[0, 0] + m[1, 1]) / 2
    si = (1 - co**2).sqrt()
    shear = -(m[0, 0] - co) / si
    height = -m[1, 0] / si
    if not height > 0:
        raise ArithmeticError("The Floquet basis must preserve orientation")
    inv = arb_mat([[1, -shear / height], [0, 1 / height]])
    all_cells = [v for panel in boxes for v in panel]
    z0, r0 = all_cells[0]
    coords = inv * arb_mat([[z0.real - 1 + arb(0, r0)], [z0.imag + arb(0, r0)]])
    u0 = acb(coords[0, 0], coords[1, 0])
    peak = arb(0)
    for i, (z, r) in enumerate(all_cells):
        j = i % len(cache)
        factor = inv * adj(core[j])
        b, forcing_tail = cache[j]
        zp, _ = local_polynomials(z, [], b, h, k)
        hp = zp - 1
        dp = 3 * zp * zp - 3
        resid = zp * zp * zp - 3 * zp + acb_poly(b)
        defect = (
            coefficient_norm(zp.derivative() + h * k * conjugate(dp) * resid)
            + h * k * coefficient_norm(dp) * forcing_tail
        )
        # The Euclidean log norm is below 10 for |z-1|<=0.01.
        err = (10 * h).exp() * (r + defect)
        if not coefficient_norm(hp) + err < arb("0.01"):
            raise ArithmeticError("Left the stationary tube")
        x = arb_poly([v.real for v in hp.coeffs()])
        y = arb_poly([v.imag for v in hp.coeffs()])
        a, bb, c, dd = arcs[j]
        u = a * x + bb * y
        v = c * x + dd * y
        X = factor[0, 0] * u + factor[0, 1] * v
        Y = factor[1, 0] * u + factor[1, 1] * v
        target = acb(0, -1) ** (i // len(cache)) * u0
        dr = X - target.real
        di = Y - target.imag
        polynomial = sum(
            abs(acb(dr[t], di[t])) for t in range(max(dr.length(), di.length()))
        )
        norm = sum(e**2 for e in factor.entries()).sqrt()
        error = norm * (
            (h * arb("55.5") * abs(k)).exp() * err + 2 * tail * coefficient_norm(hp)
        )
        peak = max(peak, ((polynomial + error) / abs(u0)).upper())
    # Residual phase difference after removing the one-period Floquet rotation.
    scale_error = (peak + abs(co.asin())).upper()
    if not scale_error < arb("0.28"):
        raise ArithmeticError(f"Four-strand tube too wide: {scale_error}")
    return {
        "template_radius": abs(u0).str(24),
        "relative_tube_width": scale_error.str(12),
        "one_turn_core_winding": core_turn.str(25),
    }


def verify(job):
    row, mesh = job
    ctx.prec = BITS
    n = PANELS
    h = 8 * arb.pi() / (n * mesh)
    k = arb(GAIN)
    cache = [drive_series(j * h, h, ORDER) for j in range(n * mesh // 4)]
    nodes = proposal(row)
    radius = arb(2) ** (-RADIUS)
    for _ in range(3):
        res, m, _, _ = follow(nodes, cache, k, h, arb(0))
        step = m.mid().solve(-res.mid())
        nodes = [
            acb((z.real + step[2 * j, 0]).mid(), (z.imag + step[2 * j + 1, 0]).mid())
            for j, z in enumerate(nodes)
        ]
    residual, m, blocks, boxes = follow(nodes, cache, k, h, radius, keep=True)
    inv = m.mid().inv()
    I = arb_mat([[int(i == j) for j in range(2 * n)] for i in range(2 * n)])
    err = I - inv * m
    q = max(
        sum(abs(err[i, j]).upper() for j in range(2 * n)).upper() for i in range(2 * n)
    )
    displacement = max(abs(x).upper() for x in inv * residual)
    if not displacement + q * radius < radius:
        raise ArithmeticError(f"Newton enclosure failed {q} {displacement} {radius}")
    M = balanced_product(blocks)
    disc = (M[0, 0] - M[1, 1]) ** 2 + 4 * M[0, 1] * M[1, 0]
    trace = M[0, 0] + M[1, 1]
    saddle = -1 + trace - M.det()
    kind = "attracting" if disc < 0 else "saddle" if saddle > 0 else "unknown"
    if kind != row["type"]:
        raise ArithmeticError(f"spectrum unresolved {kind} {disc} {saddle}")
    if not all(abs(nodes[0] - nodes[j * n // 4]) > 4 * radius for j in (1, 2, 3)):
        raise ArithmeticError("Unresolved period")
    tube = tube_bound(boxes, cache, k, h)
    return {
        "type": kind,
        **tube,
        "correction": displacement.str(12),
        "contraction": q.str(12),
        "discriminant": disc.str(12),
        "trace_minus_det_minus_one": saddle.str(12),
        "phase_zero_box_centers_relative_to_site": [
            [(nodes[j * n // 4].real - 1).str(32), nodes[j * n // 4].imag.str(32)]
            for j in range(4)
        ],
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=OUT / "satellite-boxes.json")
    p.add_argument("--substeps", type=int, default=SUBSTEPS)
    a = p.parse_args()
    if a.substeps < 1:
        p.error("Positive substeps are required")
    table = json.loads((OUT / "four-cycle.json").read_text())["records"][0]["orbits"]
    with ProcessPoolExecutor(2) as pool:
        rows = list(pool.map(verify, [(row, a.substeps) for row in table]))
    ctx.prec = BITS
    inner, outer = sorted(rows, key=lambda row: float(arb(row["template_radius"])))
    gap = (
        arb(outer["template_radius"]) * (1 - arb(outer["relative_tube_width"]))
        - arb(inner["template_radius"]) * (1 + arb(inner["relative_tube_width"]))
    ).lower()
    if not gap > 0:
        raise ArithmeticError("The two homotopies may cross")
    print("Two disjoint four-strand tubes; annular gap", gap.str(12))
    a.out.write_text(
        json.dumps(
            {
                "gain": GAIN,
                "panels": PANELS,
                "substeps": a.substeps,
                "taylor_degree": ORDER,
                "bits": BITS,
                "node_box_radius_power_of_two": -RADIUS,
                "annular_gap_lower": gap.str(16),
                "orbits": rows,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

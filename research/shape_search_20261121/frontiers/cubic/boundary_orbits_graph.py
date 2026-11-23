"""Certify a uniformly hyperbolic copy of the six-phase graph by mixed boundary data."""

import argparse
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor
from functools import cache, lru_cache
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
import itertools

import numpy as np
from cubic_material_lines import advance
from flint import acb, acb_poly, arb, arb_mat, ctx
from model import forcing
from periodic_boxes import (
    cell,
    conjugate,
    drive_series,
    local_polynomials,
)
from periodic_boxes import (
    coefficient_norm as norm,
)
from saddle_passages import KNOTS, axes, transport_cell
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[4]
SHOTS, MICRO, TILES = 32, 3, 16
BITS, DEGREE = 80, 8
EDGES = [
    (phase, site, target)
    for phase in range(3)
    for site in (-1, 1)
    for target in ((site, -site) if site == (-1) ** phase else (site,))
]


def inverse_frame(phase, site):
    u, v = axes(phase, site)
    return arb_mat([[u.real, v.real], [u.imag, v.imag]]).inv()


@cache
def cells(phase):
    first, last = KNOTS[phase : phase + 2]
    h = (last - first) * 2 * arb.pi() / (64 * SHOTS * MICRO)
    return h, [
        drive_series(2 * arb.pi() * first / 64 + h * i, h, DEGREE)
        for i in range(SHOTS * MICRO)
    ]


@lru_cache(maxsize=1)
def numerical_drive():
    source = ROOT / "results/shape_search_20261121/frontiers/cubic/certificates.json"
    row = next(
        r for r in json.loads(source.read_text())["rows"] if r["case"] == 48000001
    )
    return forcing(np.linspace(0, 2 * np.pi, 2 * 8192 + 1), row)


def predictor(phase, site, target, x, y):
    # Numerical proposals only. Each returned binary64 coordinate is treated as exact.
    drive = numerical_drive()
    start, stop = [k * 8192 // 64 for k in KNOTS[phase : phase + 2]]
    u, v = [complex(float(w.real), float(w.imag)) for w in axes(phase, site)]
    inverse = inverse_frame(phase + 1, target)
    longitudinal = np.array([float(inverse[0, i]) for i in range(2)])

    def residual(a):
        end = (
            advance(
                np.array([site + a * u + y * v]),
                drive,
                0.34,
                start,
                stop,
                2 * np.pi / 8192,
            )[0]
            - target
        )
        return longitudinal @ [end.real, end.imag] - x

    a = brentq(residual, -0.14, 0.14, xtol=1.5e-15)
    z = site + a * u + y * v
    nodes = [z]
    ticks = np.linspace(start, stop, SHOTS + 1, dtype=int)
    for first, last in itertools.pairwise(ticks):
        z = advance(np.array([z]), drive, 0.34, first, last, 2 * np.pi / 8192)[0]
        nodes.append(z)
    return np.array(nodes)


def curvature(descriptions, radius, h, gain):
    # D²Phi solves the first-variation equation with forcing D²f[DPhi, DPhi].
    # Logarithmic growth bounds both its propagator and the two forcing factors.
    second, first, dist = arb(0), arb(1), radius
    for span, az, dp, p0, defect, rounding in descriptions:
        cap = max(arb(1) / 512, 2 * dist)
        t = span + cap
        size = az + t
        dp_change = 6 * az * t + 3 * t * t
        p_max = p0 + dp * t + 3 * az * t * t + t**3 + h * 121 / 4
        mu = gain * (6 * size * p_max - max(arb(0), (dp - dp_change).lower()) ** 2)
        growth = (h * mu).exp().upper()
        peak = max(arb(1), growth)
        if not peak * (dist + defect) < cap:
            raise ArithmeticError("Second-derivative tube not enclosed")
        hessian = gain * (18 * size * (dp + dp_change) + 6 * p_max)
        second = (growth * second + h * hessian * peak**3 * first**2).upper()
        first = (growth * first).upper()
        dist = (growth * dist + peak * defect + rounding).upper()
    return second


def linearize(edge, nodes, x, y, h, gain):
    phase, site, target = edge
    n = 2 * (SHOTS + 1)
    q = arb_mat(nodes.view(float).reshape(-1, 1).tolist())
    first, last = inverse_frame(phase, site), inverse_frame(phase + 1, target)
    jac, residual, observation = arb_mat(n, n), arb_mat(n, 1), arb_mat(2, n)
    for i in range(2):
        jac[0, i], jac[n - 1, n - 2 + i] = first[1, i], last[0, i]
        observation[0, i], observation[1, n - 2 + i] = first[0, i], last[1, i]
    residual[0, 0] = jac[0, 0] * (q[0, 0] - site) + jac[0, 1] * q[1, 0] - y
    residual[n - 1, 0] = (
        jac[n - 1, n - 2] * (q[n - 2, 0] - target) + jac[n - 1, n - 1] * q[n - 1, 0] - x
    )
    base_output = observation * q - arb_mat(
        [[first[0, 0] * site], [last[1, 0] * target]]
    )
    _, forcings = cells(phase)
    descriptions = []
    for j in range(SHOTS):
        k = 2 * j
        z = acb(float(nodes[j].real), float(nodes[j].imag))
        columns = [acb(1), acb(0, 1)]
        errors = (arb(0), arb(0), [arb(0), arb(0)])
        section = []
        for b in forcings[j * MICRO : (j + 1) * MICRO]:
            path, _ = local_polynomials(z, [], b[0], h, gain)
            dp, p = 3 * path**2 - 3, path**3 - 3 * path + acb_poly(b[0])
            defect = (
                norm(path.derivative() + h * gain * conjugate(dp) * p)
                + h * gain * norm(dp) * b[1]
            )
            endpoint = path(1)
            section.append(
                (
                    norm(path - acb_poly([z])),
                    abs(z),
                    abs(3 * z * z - 3),
                    abs(z**3 - 3 * z + b[0][0]),
                    defect,
                    abs(endpoint - endpoint.mid()),
                )
            )
            z, columns, errors = cell(z, columns, errors, b, h, gain)
        descriptions.append(section)
        residual[k + 1, 0] = q[k + 2, 0] - z.real + arb(0, errors[0])
        residual[k + 2, 0] = q[k + 3, 0] - z.imag + arb(0, errors[0])
        for c, (v, e) in enumerate(zip(columns, errors[2])):
            jac[k + 1, k + c], jac[k + 2, k + c] = (
                -v.real + arb(0, e),
                -v.imag + arb(0, e),
            )
            jac[k + 1 + c, k + 2 + c] = 1
    return jac, residual, observation, base_output, descriptions


def tile(job):
    ctx.prec = BITS
    edge, i, j = job
    width, gain = arb(7) / 50, arb(17) / 50
    x = width * (arb(2 * i + 1) / TILES - 1)
    y = width * (arb(2 * j + 1) / TILES - 1)
    delta = 2 * width / TILES  # Halo reaches the centers of adjacent tiles.
    nodes = predictor(*edge, float(x), float(y))
    h, _ = cells(edge[0])
    jac, residual, observe, out0, descriptions = linearize(edge, nodes, x, y, h, gain)
    n = jac.nrows()
    midpoint = np.array([[float(jac[a, b].mid()) for b in range(n)] for a in range(n)])
    B = arb_mat(np.linalg.inv(midpoint).tolist())
    identity = arb_mat(np.eye(n, dtype=int).tolist())
    # Mixed parameters are (target x, source y), in that order.
    S = arb_mat([[B[k, n - 1], B[k, 0]] for k in range(n)])
    offset = B * residual
    displacement = [
        (abs(S[k, 0]) + abs(S[k, 1])) * delta + abs(offset[k, 0]) for k in range(n)
    ]
    radii = [(2 * d + arb(1) / 10**6).upper() for d in displacement]
    for _ in range(4):
        box_jac = arb_mat(jac)
        for k, desc in enumerate(descriptions):
            radius = (radii[2 * k] ** 2 + radii[2 * k + 1] ** 2).sqrt().upper()
            error = arb(0, (curvature(desc, radius, h, gain) * radius).upper())
            for a in (2 * k + 1, 2 * k + 2):
                for b in (2 * k, 2 * k + 1):
                    box_jac[a, b] += error
        E = identity - B * box_jac
        image = [
            sum(abs(E[k, b]) * radii[b] for b in range(n)).upper() for k in range(n)
        ]
        needed = [(a + b).upper() for a, b in zip(displacement, image)]
        if all(a < b for a, b in zip(needed, radii)):
            break
        radii = [max(r, (arb(11) * a / 10).upper()) for r, a in zip(radii, needed)]
    else:
        raise ArithmeticError(f"No uniform Krawczyk inclusion at {job}")
    q = max((a / r).upper() for a, r in zip(image, radii))
    center_error = (
        max((abs(offset[k, 0]) / radii[k]).upper() for k in range(n)) / (1 - q)
    ).upper()
    # The inclusion gives |S| 1 + |E| radii/delta < radii/delta. Thus each
    # row norm of Dz/Dp is <= radii/delta in Dz/Dp = S + E Dz/Dp.
    constant, remainder = observe * S, observe * E
    derivative = [
        (
            abs(constant[k, 0])
            + abs(constant[k, 1])
            + sum(abs(remainder[k, b]) * radii[b] for b in range(n)) / delta
        ).upper()
        for k in range(2)
    ]
    if not max(derivative) < 1:
        raise ArithmeticError(f"Cross map lacks contraction at {job}")
    for k in range(2):
        value = abs(out0[k, 0]) + sum(
            abs(observe[k, b]) * radii[b] * center_error for b in range(n)
        )
        if not value + derivative[k] * delta < width:
            raise ArithmeticError(f"Cross map did not return to the square at {job}")
    # Only dyadic enclosures cross the worker boundary. Check their rounding directions.
    small = [
        float(np.nextafter(float((r * center_error).upper()), np.inf)) for r in radii
    ]
    large = [float(np.nextafter(float(r.lower()), -np.inf)) for r in radii]
    bound = float(np.nextafter(float(max(derivative)), np.inf))
    if not (
        arb(bound) > max(derivative)
        and all(arb(a) > r * center_error for a, r in zip(small, radii))
        and all(arb(a) < r for a, r in zip(large, radii))
    ):
        raise ArithmeticError("Enclosure could not be exported")
    return edge, i, j, (nodes.view(float).tolist(), small, large, bound)


def fits(first, second):
    # At first's parameter center the root lies in second's uniqueness box.
    a, rad, _, _ = first
    b, _, box, _ = second
    return all(
        abs(arb(x) - arb(y)) + arb(r) < arb(R) for x, y, r, R in zip(a, b, rad, box)
    )


def glue(edge, tiles):
    seams = 0
    for i in range(TILES):
        for j in range(TILES):
            for ni, nj in ((i + 1, j), (i, j + 1)):
                if ni == TILES or nj == TILES:
                    continue
                a, b = tiles[i, j], tiles[ni, nj]
                if not (fits(a, b) or fits(b, a)):
                    raise ArithmeticError(
                        f"Adjacent branches not identified: {edge, i, j, ni, nj}"
                    )
                seams += 1
    if edge[1] == edge[2]:
        fixed = [v for _ in range(SHOTS + 1) for v in (edge[1], 0)]
        center_tiles = [
            tiles[i, j]
            for i in (TILES // 2 - 1, TILES // 2)
            for j in (TILES // 2 - 1, TILES // 2)
        ]
        if not any(
            all(abs(arb(x) - arb(y)) < arb(r) for x, y, r in zip(fixed, v[0], v[2]))
            for v in center_tiles
        ):
            raise ArithmeticError("Stationary orbit did not match its selected branch")
    return {
        "phase": edge[0],
        "source": edge[1],
        "target": edge[2],
        "tiles": len(tiles),
        "seams": seams,
        "derivative_inf_norm_upper": max(v[3] for v in tiles.values()),
    }


def match_archived_saddles(tiles):
    directory = ROOT / "results/cubic-forced-orbits"
    proof = json.loads((directory / "orbits.json").read_text())
    path = directory / "shooting-nodes.npz"
    if (
        hashlib.sha256(path.read_bytes()).hexdigest() != proof["node_file_sha256"]
        or proof["box_half_width"] != "2^-40"
        or proof["phase_nodes_per_turn"] != 64
    ):
        raise ArithmeticError(
            "Shooting archive differs from its periodic-orbit certificate"
        )
    nodes = np.load(path)["nodes"]
    step = 2 * arb.pi() / (64 * SHOTS)
    forcing_data = [drive_series(i * step, step, DEGREE) for i in range(64 * SHOTS)]
    width, gain = arb(7) / 50, arb(17) / 50
    core = width / TILES
    centers = [width * (arb(2 * i + 1) / TILES - 1) for i in range(TILES)]

    def coordinate(disk, phase, site, axis):
        z, radius = disk
        A = inverse_frame(phase, site)
        return (
            A[axis, 0] * (z.real - site)
            + A[axis, 1] * z.imag
            + arb(0, radius * (A[axis, 0] ** 2 + A[axis, 1] ** 2).sqrt())
        )

    matched = []
    for row in proof["orbits"]:
        if "saddle" not in row["type"]:
            continue
        orbit = row["orbit"]
        disks = []
        for j, start in enumerate(nodes[orbit]):
            z = acb(float(start.real), float(start.imag))
            radius = arb(2).sqrt() * arb(2) ** -40
            for f in forcing_data[j * SHOTS : (j + 1) * SHOTS]:
                disks.append((z, radius))
                z, radius = transport_cell(z, radius, f, step, gain)
        for phase, (first, last) in enumerate(itertools.pairwise(KNOTS)):
            samples = [
                disks[(SHOTS * first + i * (last - first)) % len(disks)]
                for i in range(SHOTS + 1)
            ]
            if phase >= 3:
                samples = [(-z.conjugate(), r) for z, r in samples]
            site, target = [1 if samples[i][0].real > 0 else -1 for i in (0, -1)]
            p = phase % 3
            x, y = (
                coordinate(samples[-1], p + 1, target, 0),
                coordinate(samples[0], p, site, 1),
            )
            if not (abs(x) < width and abs(y) < width):
                raise ArithmeticError("Archived saddle outside mixed-data square")
            hit = False
            for i, xc in enumerate(centers):
                for j, yc in enumerate(centers):
                    if abs(x - xc) > core or abs(y - yc) > core:
                        continue
                    if not (abs(x - xc) < 2 * core and abs(y - yc) < 2 * core):
                        raise ArithmeticError("Unresolved parameter tile")
                    node, _, radius, _ = tiles[p, site, target][i, j]
                    for k, (z, r) in enumerate(samples):
                        for c, value in enumerate((z.real, z.imag)):
                            if not abs(value - arb(node[2 * k + c])) + r < arb(
                                radius[2 * k + c]
                            ):
                                raise ArithmeticError(
                                    f"Saddle {orbit} fell outside its uniqueness box"
                                )
                    hit = True
            if not hit:
                raise ArithmeticError("No matching mixed-data tile")
        matched.append(orbit)
    if len(matched) != 18:
        raise ArithmeticError("Expected 18 archived saddle orbits")
    return matched


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results/cubic-saddle-passages/hyperbolicity.json",
    )
    args = parser.parse_args()
    ctx.prec = BITS
    jobs = [(edge, i, j) for edge in EDGES for i in range(TILES) for j in range(TILES)]
    tiles = {edge: {} for edge in EDGES}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for edge, i, j, result in pool.map(tile, jobs, chunksize=4):
            tiles[edge][i, j] = result
    records = [glue(edge, tiles[edge]) for edge in EDGES]
    matched = match_archived_saddles(tiles)
    packet = {
        "case": 48000001,
        "gain": "17/50",
        "half_width": "7/50",
        "phase_knots_over_64": KNOTS,
        "shots": SHOTS,
        "cells_per_shot": MICRO,
        "taylor_degree": DEGREE,
        "arb_bits": BITS,
        "tiles_per_axis": TILES,
        "branches": records,
        "archived_saddles_in_graph": matched,
    }
    args.out.write_text(json.dumps(packet, indent=2) + "\n", encoding="utf-8")
    print(
        "Certified mixed derivative bound:",
        max(r["derivative_inf_norm_upper"] for r in records),
    )

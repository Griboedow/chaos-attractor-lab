"""Certify three basin exits from each moving one-period saddle of the driven cubic."""

import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
from flint import acb, acb_poly, arb, arb_mat, ctx
from periodic_boxes import cell, conjugate, drive_series, local_polynomials
from periodic_boxes import coefficient_norm as norm
from saddle_passages import transport_cell

ROOT = Path(__file__).resolve().parents[4]
PERIODIC = ROOT / "results/cubic-forced-orbits"
CUBIC = ROOT / "results/shape_search_20261121/frontiers/cubic"
SUBSTEPS = 32
DEGREE = 10
BITS = 256
# Phase in 64ths, then (signed displacement, full turns) for sinks A, B and C.
# The remaining fraction of a turn always returns the test to phase zero.
STARTS = {
    5: (0, ("3.73e-17", 2), ("-1.02e-17", 2), ("2.1e-16", 2)),
    6: (0, ("9.58e-31", 2), ("-2.62e-31", 2), ("5.39e-30", 2)),
    7: (32, ("-9.58e-31", 2), ("2.62e-31", 2), ("-5.39e-30", 2)),
    8: (0, ("3.71e-24", 2), ("-1.02e-24", 2), ("2.09e-23", 2)),
    9: (53, ("-4.28e-23", 2), ("-4.28e-24", 2), ("1.02e-24", 2)),
    10: (0, ("7.5e-32", 2), ("4.22e-35", 3), ("-1.78e-31", 2)),
    11: (53, ("-2.32e-30", 2), ("-2.32e-31", 2), ("5.51e-32", 2)),
    12: (24, ("1.31e-29", 2), ("2.32e-31", 2), ("-3.58e-31", 2)),
    13: (0, ("5.72e-30", 2), ("-1.02e-30", 2), ("-4.29e-32", 3)),
    14: (60, ("3.22e-30", 3), ("-2.42e-31", 2), ("5.73e-35", 3)),
    15: (44, ("-1.25e-24", 2), ("2.96e-23", 2), ("3.42e-24", 2)),
    16: (0, ("2.22e-24", 2), ("-6.08e-25", 2), ("1.25e-23", 2)),
    17: (32, ("-1.17e-23", 2), ("3.19e-24", 2), ("-6.56e-23", 2)),
    18: (0, ("1.17e-23", 2), ("-3.19e-24", 2), ("6.56e-23", 2)),
    19: (44, ("-3.27e-31", 2), ("6.71e-30", 2), ("8.95e-31", 2)),
    20: (12, ("3.27e-31", 2), ("-6.71e-30", 2), ("-8.95e-31", 2)),
}


def field_bound(path, b0, error, h, gain):
    center = path[0]
    span = norm(path - acb_poly([center])) + error
    size = abs(center) + span
    dp = abs(3 * center**2 - 3)
    dp_min = max(arb(0), (dp - 6 * abs(center) * span - 3 * span**2).lower())
    p_max = (
        abs(center**3 - 3 * center + b0)
        + dp * span
        + 3 * abs(center) * span**2
        + span**3
        + h * 121 / 4
    )
    return (gain * (-(dp_min**2) + 6 * size * p_max)).upper(), size


def column_matrix(columns, errors):
    matrix = arb_mat(2, 2)
    for j, (v, error) in enumerate(zip(columns, errors)):
        matrix[0, j] = v.real + arb(0, error)
        matrix[1, j] = v.imag + arb(0, error)
    return matrix


def reference_flow(nodes, proof, phase, h, gain, substeps, degree):
    initial_radius = (
        arb(2).sqrt()
        * arb(proof["center_correction_bound"]).upper()
        / (1 - arb(proof["contraction_bound"]).upper())
    )
    if not initial_radius < arb(2).sqrt() * arb(2) ** -40:
        raise ArithmeticError("Unexpected interval Newton bound")
    matrix = arb_mat([[1, 0], [0, 1]])
    lipschitz, quadratic = arb(1), arb(0)
    tube = arb(1) / 1024
    cells = []
    for j in range(64):
        k = (phase + j) % 64
        z = acb(float(nodes[k].real), float(nodes[k].imag))
        columns = [acb(1), acb(0, 1)]
        errors = (arb(0), initial_radius, [arb(0), arb(0)])
        for i in range(substeps):
            b = drive_series((k * substeps + i) * h, h, degree)
            path, _ = local_polynomials(z, [], b[0], h, gain)
            dp, p = 3 * path**2 - 3, path**3 - 3 * path + acb_poly(b[0])
            defect = (
                norm(path.derivative() + h * gain * conjugate(dp) * p)
                + h * gain * norm(dp) * b[1]
            )
            mu, _ = field_bound(path, b[0][0], tube, h, gain)
            error = (max(arb(1), (h * mu).exp()) * (errors[1] + defect)).upper()
            if not error < tube:
                raise ArithmeticError("Reference orbit left its Taylor tube")
            mu, size = field_bound(path, b[0][0], tube + error, h, gain)
            if not size < 4:
                raise ArithmeticError("Hessian bound left its domain")
            growth = (h * mu).exp().upper()
            peak = max(arb(1), growth)
            derivative = column_matrix(columns, errors[2]) * matrix
            frobenius_squared = sum(abs(entry) ** 2 for entry in derivative)
            # Hessian norm < 2000 on |z| <= 4. Integrate the Taylor remainder
            # along the actual variational flow, retaining contraction between cells.
            quadratic = (
                growth * quadratic + h * peak**3 * 1000 * frobenius_squared
            ).upper()
            lipschitz = (lipschitz * peak).upper()
            cells.append((b, path, error, dp, p))
            z, columns, errors = cell(z, columns, errors, b, h, gain)
        matrix = column_matrix(columns, errors[2]) * matrix
    trace = matrix[0, 0] + matrix[1, 1]
    if not abs(trace) > 2:
        raise ArithmeticError("Unresolved hyperbolic saddle")
    # Divergence -2*kappa*|p'|^2 gives 0 < det <= 1.
    sign = 1 if trace > 0 else -1
    eigenvalue = (trace + sign * (trace**2 - 4 * arb(0.5, 0.5)).sqrt()) / 2
    vector = acb(matrix[0, 1], eigenvalue - matrix[0, 0])
    vector /= abs(vector)
    if vector.real < 0:
        vector = -vector
    if not (vector.real > 0 and eigenvalue**2 > lipschitz):
        raise ArithmeticError("Unresolved unstable direction")
    error_factor = (quadratic / (eigenvalue**2 - lipschitz)).upper()
    return cells, vector, eigenvalue, lipschitz, error_factor


def difference_cell(w, radius, reference, h, gain):
    b, z, error, dp, p = reference
    total, _ = local_polynomials(z[0] + w, [], b[0], h, gain)
    v = total - z
    size = norm(v)
    change = 3 * v * (2 * z + v)
    residual_change = v * (dp + v * (3 * z + v))
    rhs = conjugate(dp + change) * residual_change + conjugate(change) * p
    defect = norm(v.derivative() + h * gain * rhs) + h * gain * norm(change) * b[1]
    # The derivative of f(z+v)-f(z) with respect to z is bounded by 2000|v|.
    defect += h * 2000 * error * size
    cap = max(arb(1) / 10**8, 2 * radius)
    mu, disk_size = field_bound(total, b[0][0], error + cap, h, gain)
    growth = (h * mu).exp().upper()
    peak = max(arb(1), growth)
    if not (
        disk_size < 4
        and norm(z) + error + size + cap < 4
        and peak * (radius + defect) < cap
    ):
        raise ArithmeticError("Perturbed trajectory left its enclosure")
    end = v(1)
    return end.mid(), (growth * radius + peak * defect + abs(end - end.mid())).upper()


def propagate(scale, steps, cells, vector, lipschitz, error_factor, h, gain):
    s = arb(scale)
    if not abs(s) + error_factor * s**2 < 1 / (1024 * lipschitz):
        raise ArithmeticError("Start outside local return-map bound")
    initial = s * vector
    z = initial.mid()
    radius = (abs(initial - z) + error_factor * s**2).upper()
    relative = True
    for index in range(steps):
        ref = cells[index % len(cells)]
        if relative and abs(z) > arb(1) / 200:
            total = z + ref[1][0]
            radius = (radius + ref[2] + abs(total - total.mid())).upper()
            z, relative = total.mid(), False
        if relative:
            z, radius = difference_cell(z, radius, ref, h, gain)
        else:
            enclosed = transport_cell(z, radius, ref[0], h, gain)
            if enclosed is None:
                raise ArithmeticError(
                    "Exit disk left the absolute-coordinate enclosure"
                )
            z, radius = enclosed
    if relative:
        ref = cells[steps % len(cells)]
        z, radius = z + ref[1][0], radius + ref[2]
    return z, radius


def verify(orbit, substeps=SUBSTEPS, degree=DEGREE):
    ctx.prec = BITS
    gain = arb(17) / 50
    if not gain * (18 * 4 * 51 + 6 * (4**3 + 12 + arb(29) / 4)) < 2000:
        raise ArithmeticError("Hessian bound failed")
    h = 2 * arb.pi() / (64 * substeps)
    info = json.loads((PERIODIC / "orbits.json").read_text())
    path = PERIODIC / "shooting-nodes.npz"
    if (
        hashlib.sha256(path.read_bytes()).hexdigest() != info["node_file_sha256"]
        or info["gain"] != [17, 50]
        or info["case"] != 48000001
        or info["phase_nodes_per_turn"] != 64
        or info["box_half_width"] != "2^-40"
    ):
        raise ValueError("Unexpected periodic-orbit data")
    proof = next(q for q in info["orbits"] if q["orbit"] == orbit)
    with np.load(path) as source:
        nodes = source["nodes"][orbit]
    survey = next(
        q
        for q in json.loads((CUBIC / "gradient.json").read_text())
        if q["case"] == info["case"]
    )
    targets = next(q for q in survey["runs"] if q["gain"] == 0.34)["starts_last_period"]
    checks = json.loads((CUBIC / "finite-gains.json").read_text())["records"]
    trap = next(q for q in checks if q["gain"] == info["gain"])["checks"][-1]
    if not trap["success"]:
        raise ValueError("Missing trapping certificate")
    phase, *seeds = STARTS[orbit]
    cells, vector, eigenvalue, lipschitz, factor = reference_flow(
        nodes, proof, phase, h, gain, substeps, degree
    )
    exits = []
    for target, (scale, turns) in enumerate(seeds):
        steps = (64 * turns + (-phase) % 64) * substeps
        z, radius = propagate(scale, steps, cells, vector, lipschitz, factor, h, gain)
        distance = (abs(z - acb(*targets[target])) + radius).upper()
        if not distance < arb(float(trap["start_radius"])):
            raise ArithmeticError(f"Orbit {orbit} missed sink {target}: {distance}")
        exits.append(
            {
                "target": target,
                "scale": scale,
                "elapsed_over_64": steps // substeps,
                "distance_upper": distance.str(8),
            }
        )
    print(f"Orbit {orbit:02d}: three basin exits", flush=True)
    return {
        "orbit": orbit,
        "phase_over_64": phase,
        "unstable_multiplier": eigenvalue.str(8),
        "log_lipschitz_upper": lipschitz.log().str(8),
        "quadratic_factor_upper": factor.str(8),
        "exits": exits,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--substeps", type=int, default=SUBSTEPS)
    parser.add_argument("--degree", type=int, default=DEGREE)
    parser.add_argument(
        "--orbit", choices=STARTS, type=int, nargs="+", default=sorted(STARTS)
    )
    default_output = ROOT / "results/cubic-saddle-passages/boundary-orbits.json"
    parser.add_argument("--out", type=Path, default=default_output)
    args = parser.parse_args()
    if args.substeps < 1 or args.degree < 2 or args.workers < 1:
        parser.error(
            "Positive substep/worker counts and Taylor degree at least two are required"
        )
    if args.out.resolve() == default_output and set(args.orbit) != set(STARTS):
        parser.error("Use --out to save a partial orbit check outside the full packet")
    run = partial(verify, substeps=args.substeps, degree=args.degree)
    if args.workers == 1:
        records = list(map(run, args.orbit))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            records = list(pool.map(run, args.orbit))
    receipt = {
        "case": 48000001,
        "gain": "17/50",
        "bits": BITS,
        "degree": args.degree,
        "cells_per_period": 64 * args.substeps,
        "moving_saddles": records,
    }
    args.out.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")

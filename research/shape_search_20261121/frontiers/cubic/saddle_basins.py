"""Enclose three forward images of the local unstable manifold of z=1."""

import argparse
import hashlib
import json
from pathlib import Path

from flint import acb, acb_poly, arb, ctx
from periodic_boxes import cell, coefficient_norm, conjugate, drive_series

ROOT = Path(__file__).resolve().parents[4]
DATA = ROOT / "results/shape_search_20261121/frontiers/cubic"


def linear_cell(start, a, h, tail, gain):
    """Taylor path and integrated defect for v' = a(t) conj(v)."""
    v = [start]
    for k in range(len(a) - 1):
        product = sum((a[j] * v[k - j].conjugate() for j in range(k + 1)), acb(0))
        v.append(h * product / (k + 1))
    vp = acb_poly(v)
    defect = coefficient_norm(vp.derivative() - h * acb_poly(a) * conjugate(vp))
    defect += h * 6 * gain * tail * coefficient_norm(vp)
    return vp, defect


def eigenvector(forcings, h, gain, reverse=False):
    # At z=1: -6*k*(b-2) forward, -6*k*(b+2) in reversed time since b is odd.
    columns = [acb(1), acb(0, 1)]
    errors = [arb(0), arb(0)]
    for b, tail in forcings:
        offset = -2 if reverse else 2
        a = [-6 * gain * (x - (offset if i == 0 else 0)) for i, x in enumerate(b)]
        peak = (h * (abs(a[0]) + 6 * gain * h * 121 / 4)).exp().upper()
        for c, start in enumerate(columns):
            vp, defect = linear_cell(start, a, h, tail, gain)
            end = vp(1)
            center = end.mid()
            errors[c] = (peak * (errors[c] + defect) + abs(end - center)).upper()
            columns[c] = center
    c0, c1 = [v + acb(arb(0, e), arb(0, e)) for v, e in zip(columns, errors)]
    trace = c0.real + c1.imag
    # det M = 1 because div F vanishes along this exact orbit.
    eigenvalue = (trace - (trace**2 - 4).sqrt()) / 2
    assert eigenvalue < -(arb(2) ** 63)
    u = acb(c1.real, eigenvalue - c0.real)
    return eigenvalue, u / abs(u)


def local_estimates(gain, reverse=False):
    delta = arb(1) / 64
    # RMS of b-2 uses orthogonality of harmonics. Gronwall on |z-1|<=delta.
    rms = (arb(529) / 32).sqrt()
    h = 2 * arb.pi() * 6 * gain * (1 + delta) * (rms + 3 * delta**2 + delta**3)
    size = 1 + delta
    dp = 6 * delta + 3 * delta**2
    if reverse:
        h += 2 * arb.pi() * gain * dp**2
    assert h.exp() < arb(2) ** 78
    pmax = arb(37) / 4 + 3 * delta**2 + delta**3
    jacobian_lipschitz = gain * (18 * size * dp + 6 * pmax)
    assert arb.pi() * jacobian_lipschitz < 256
    quadratic = arb(2) ** 164
    return quadratic / (arb(2) ** 126 - arb(2) ** 78)


def half_turn(sinks, trapping):
    """Locate the three reflected certified sink boxes in the phase-zero trapping disks."""
    import numpy as np

    folder = ROOT / "results/cubic-forced-orbits"
    proof = json.loads((folder / "orbits.json").read_text())
    path = folder / "shooting-nodes.npz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == proof["node_file_sha256"]
    assert proof["box_half_width"] == "2^-40"
    nodes = np.load(path)["nodes"]
    radius = arb(2).sqrt() * arb(2) ** -40
    permutation = []
    for orbit in proof["orbits"]:
        if orbit["type"] != "sink":
            continue
        row = nodes[orbit["orbit"]]
        middle = row[len(row) // 2]
        zero = acb(float(row[0].real), float(row[0].imag))
        assert abs(zero - sinks[len(permutation)]) + radius < trapping
        reflected = acb(-float(middle.real), float(middle.imag))
        hits = [
            i
            for i, sink in enumerate(sinks)
            if abs(reflected - sink) + radius < trapping
        ]
        assert len(hits) == 1
        permutation.append(hits[0])
    assert sorted(permutation) == [0, 1, 2]
    return permutation


def main(args):
    ctx.prec = 512
    gain = arb(17) / 50
    h = 2 * arb.pi() / args.cells
    forcing = [drive_series(j * h, h, args.degree) for j in range(args.cells)]
    eigenvalue, unit = eigenvector(forcing, h, gain)
    print(
        "eigen",
        eigenvalue.str(12),
        unit.str(12),
        "direction radius",
        abs(unit - unit.mid()).upper().str(5),
        flush=True,
    )
    coefficient = local_estimates(gain)
    packet = json.loads((DATA / "finite-gains.json").read_text())
    check = next(row for row in packet["records"] if row["gain"] == [17, 50])["checks"][
        -1
    ]
    assert check["success"]
    trapping = arb(float(check["start_radius"]))
    survey = next(
        x
        for x in json.loads((DATA / "gradient.json").read_text())
        if x["case"] == 48000001
    )
    sink = next(row for row in survey["runs"] if row["gain"] == 0.34)[
        "starts_last_period"
    ]
    targets = [acb(x, y) for x, y in sink]
    records = []
    for scale in ("1e-73", "-1e-73", "1e-70"):
        s = arb(scale)
        assert abs(s) + coefficient * s * s < arb(2) ** -84
        z = 1 + s * unit
        center = z.mid()
        error = (abs(z - center) + coefficient * s * s).upper()
        initial_error = error
        for turn in range(4):
            for forcing_cell in forcing:
                center, _, errors = cell(
                    center, [], (error, error, []), forcing_cell, h, gain
                )
                error = errors[1]
            print(scale, turn + 1, center.str(10), error.str(5), flush=True)
        distance = [(abs(center - target) + error).upper() for target in targets]
        hits = [i for i, x in enumerate(distance) if x < trapping]
        assert len(hits) == 1
        records.append(
            {
                "scale": scale,
                "target": hits[0],
                "turns": 4,
                "initial_disk_radius": str(initial_error),
                "final_disk_radius": str(error),
                "target_distance_upper": str(distance[hits[0]]),
                "endpoint": str(center),
            }
        )
    assert len({r["target"] for r in records}) == 3
    result = {
        "case": 48000001,
        "gain": [17, 50],
        "saddle": [1, 0],
        "precision_bits": ctx.prec,
        "cells_per_turn": args.cells,
        "degree": args.degree,
        "eigenvalue": str(eigenvalue),
        "unit_unstable_direction": str(unit),
        "trapping_radius": check["start_radius"],
        "connections": records,
        "half_turn_sink_permutation": half_turn(targets, trapping),
        "unit_midpoint": [float(unit.real.mid()), float(unit.imag.mid())],
        "sink_centers": sink,
    }
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


def stable_crossing(args):
    """Transport an entire stable arc through the phase-zero passage box."""
    from saddle_passages import axes

    ctx.prec = 512
    gain = arb(17) / 50
    h = 2 * arb.pi() / args.cells
    forcing = [drive_series(j * h, h, args.degree) for j in range(args.cells)]
    eigenvalue, unit = eigenvector(forcing, h, gain, reverse=True)
    if unit.imag < 0:
        unit = -unit
    scale = arb(1) / (5 * eigenvalue)
    coefficient = local_estimates(gain, reverse=True)
    previous = abs(scale / eigenvalue)
    assert previous + coefficient * previous**2 < arb(2) ** -84
    initial = scale * unit
    center = initial.mid()
    error = (abs(initial - center) + coefficient * scale**2).upper()
    for b, tail in forcing:
        a = [-6 * gain * (bj + (2 if j == 0 else 0)) for j, bj in enumerate(b)]
        path, defect = linear_cell(center, a, h, tail, gain)
        r = coefficient_norm(path)
        drive_max = abs(b[0] + 2) + h * 121 / 4
        # For every real |t|<=1, t*path is a common nominal solution up to
        # the linear defect plus the exact nonlinear remainder at |w|<=r.
        defect += (
            h * gain * (3 * r**2 * drive_max + (6 * r + 3 * r**2) * (3 + r) * r**2)
        )
        cap = (2 * (error + defect) + arb(2) ** -256).upper()
        size = r + cap
        dp = 6 * size + 3 * size**2
        p_max = drive_max + 3 * size**2 + size**3
        growth = (h * gain * (dp**2 + 6 * (1 + size) * p_max)).exp().upper()
        peak = (growth * (error + defect)).upper()
        assert peak < cap
        end = path(1)
        center = end.mid()
        error = (peak + abs(end - center)).upper()
    horizontal, vertical = axes(0, 1)
    determinant = (horizontal.conjugate() * vertical).imag
    assert determinant > 0
    x = (center.conjugate() * vertical).imag / determinant
    y = (horizontal.conjugate() * center).imag / determinant
    x_bound = (abs(x) + error * abs(vertical) / determinant).upper()
    y_bound = (y - error * abs(horizontal) / determinant).lower()
    if not x_bound < arb(7) / 50 < y_bound:
        raise ArithmeticError(
            f"Unresolved arc: end={center.str(8)} error={error.str(6)} x={x_bound.str(6)} y={y_bound.str(6)}"
        )
    result = {
        "case": 48000001,
        "gain": [17, 50],
        "phase": 0,
        "site": 1,
        "bits": ctx.prec,
        "cells": args.cells,
        "degree": args.degree,
        "box_half_width": "7/50",
        "expanding_reverse_eigenvalue": str(eigenvalue),
        "reverse_scale": str(scale),
        "arc_nominal_endpoint": str(center),
        "uniform_euclidean_error": str(error),
        "absolute_horizontal_upper": str(x_bound),
        "positive_endpoint_vertical_lower": str(y_bound),
        "negative_endpoint_vertical_upper": str(-y_bound),
    }
    print(
        "stable arc: |x| <=",
        x_bound.str(7),
        "; signed y >=",
        y_bound.str(7),
        flush=True,
    )
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path)
    p.add_argument("--degree", type=int, default=28)
    p.add_argument(
        "--cells",
        type=int,
        help="Cells per turn (512 for exits; 2048 for the stable arc).",
    )
    p.add_argument(
        "--stable-crossing",
        action="store_true",
        help="Verify the phase-zero stable arc.",
    )
    args = p.parse_args()
    if args.cells is None:
        args.cells = 2048 if args.stable_crossing else 512
    (stable_crossing if args.stable_crossing else main)(args)

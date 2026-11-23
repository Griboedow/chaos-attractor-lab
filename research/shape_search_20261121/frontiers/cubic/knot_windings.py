"""Enclose the meridional winding of each resonant stationary variational orbit."""

import argparse
import json
from pathlib import Path

from elliptic_sites import OUT, eye, local_flow
from flint import arb, arb_mat, ctx


def winding(row, cells, degree):
    gain = arb(row["gain_center"]) + arb(0, arb(row["gain_radius"]).upper())
    step = 2 * arb.pi() / cells
    complex_radius = 8 * step
    bound = abs(gain) * (
        12
        + arb(3) / 2 * complex_radius.cosh()
        + 24 * (3 * complex_radius).cosh()
        + 18 * (6 * complex_radius).cosh()
    )
    tail = (complex_radius * bound).exp() * arb(8) ** (-degree) / 7
    matrices = [eye()]
    for j in range(cells):
        local, _ = local_flow(j * step, step, gain, 1, degree, tail, arb(0), False)
        matrices.append(local * matrices[-1])
    # For v'=c(t)*conj(v), |d arg(v)/dt| <= |c(t)|. This eliminates
    # aliasing between consecutive sample phases, independently of |v|.
    step_angle = step * abs(gain) * arb("55.5")
    if not step_angle < arb.pi() / 2:
        raise ArithmeticError("Need a finer angular mesh")
    first = arb_mat([[1], [0]])
    total = arb(0)
    for _ in range(row["turns"]):
        previous = first
        for matrix in matrices[1:]:
            new = matrix * first
            x, y = previous[0, 0], previous[1, 0]
            X, Y = new[0, 0], new[1, 0]
            dot = x * X + y * Y
            if not dot > 0:
                raise ArithmeticError(
                    "Cannot resolve the short arc between the vectors"
                )
            total += ((x * Y - y * X) / dot).atan()
            previous = new
        first = previous
    rotations = total / (2 * arb.pi())
    integer = round(float(rotations.mid()))
    if not abs(rotations - integer) < arb(1) / 8:
        raise ArithmeticError("Cannot isolate an integer winding")
    return {
        "turns": row["turns"],
        "gain_center": row["gain_center"],
        "winding": integer,
        "enclosed_rotations": rotations.str(20),
        "max_angle_per_cell": step_angle.upper().str(10),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=96)
    parser.add_argument("--degree", type=int, default=60)
    parser.add_argument("--out", type=Path, default=OUT / "knot-windings.json")
    args = parser.parse_args()
    ctx.prec = 256
    data = json.loads((OUT / "resonances.json").read_text(encoding="utf-8"))
    rows = [winding(row, args.cells, args.degree) for row in data["resonances"]]
    print([(x["turns"], x["winding"]) for x in rows])
    args.out.write_text(
        json.dumps(
            {
                "site": 1,
                "bits": 256,
                "cells_per_turn": args.cells,
                "degree": args.degree,
                "windings": rows,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

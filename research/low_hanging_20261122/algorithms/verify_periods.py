"""Verify separation and least periods for the archived R1051 covering windows."""

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from flint import arb, ctx, fmpq

ROOT = Path(__file__).resolve().parents[3]
INPUT = ROOT / "results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json"
RECEIPT = ROOT / "results/low_hanging_20261122/algorithms/certificate.json"


def ball(q):
    return arb(fmpq(q.numerator, q.denominator))


def dot(v, w):
    return sum(a * b for a, b in zip(v, w))


def cross(a, b):
    return tuple(
        a[(j + 1) % 3] * b[(j + 2) % 3] - a[(j + 2) % 3] * b[(j + 1) % 3] for j in range(3)
    )


def separated_windows(centers, columns):
    """Exact projections onto face normals and cross products of edge directions."""
    margins = []
    for i, j in combinations(range(7), 2):
        delta = [a - b for a, b in zip(centers[i], centers[j])]
        edges = [*columns[i], *columns[j]]
        best = Fraction(0)
        for left, right in combinations(edges, 2):
            axis = cross(left, right)
            length = sum(map(abs, axis))
            if length:
                gap = abs(dot(axis, delta)) - sum(abs(dot(axis, edge)) for edge in edges)
                best = max(best, gap / length)
        assert best > Fraction(1, 2000), (i, j)
        margins.append({"pair": [i, j], "linf_gap_lower_millionths": int(10**6 * best)})
    return margins


def step(x, params, tangents=()):
    beta, r, cs, sn, center = params
    dx, dy, dz = [xi - ci for xi, ci in zip(x, center)]
    v = dy * cs + dz * sn
    norm_b = (dx**2 + v**2).sqrt()
    if not norm_b > 0:
        return None
    bx, by, bz = (
        center[0] + r * dx / norm_b,
        center[1] + r * cs * v / norm_b,
        center[2] + r * sn * v / norm_b,
    )
    rx, ry = 2 * bx - x[0], 2 * by - x[1]
    norm_a = (rx**2 + ry**2).sqrt()
    if not norm_a > 0:
        return None
    y = [
        beta * x[0] + beta * rx / norm_a + (1 - 2 * beta) * bx,
        beta * x[1] + beta * ry / norm_a + (1 - 2 * beta) * by,
        beta * x[2] + (1 - 2 * beta) * bz,
    ]
    new_tangents = []
    for ux, uy, uz in tangents:
        # Each circle projection differentiates along its in-plane tangent.
        scale_b = r * (-v * ux + dx * (cs * uy + sn * uz)) / norm_b**3
        dbx, dby, dbz = -v * scale_b, dx * cs * scale_b, dx * sn * scale_b
        scale_a = (-ry * (2 * dbx - ux) + rx * (2 * dby - uy)) / norm_a**3
        new_tangents.append(
            [
                beta * (ux - ry * scale_a) + (1 - 2 * beta) * dbx,
                beta * (uy + rx * scale_a) + (1 - 2 * beta) * dby,
                beta * uz + (1 - 2 * beta) * dbz,
            ]
        )
    return y, new_tangents


def enclosed_box(center, cols, local):
    coordinates = [arb(ball(mid), ball(radius)) for mid, radius in local]
    return [
        ball(c) + sum(ball(edge[k]) * t for edge, t in zip(cols, coordinates))
        for k, c in enumerate(center)
    ]


def phase_separation(centers, columns, data):
    beta, r, theta, *center = [ball(Fraction(v)) for v in data["params"]]
    params = beta, r, theta.cos(), theta.sin(), center
    full = [(Fraction(0), Fraction(1)) for _ in range(3)]
    windows = [enclosed_box(c, e, full) for c, e in zip(centers, columns)]
    hull = [windows[0][i] for i in range(3)]
    for box in windows[1:]:
        hull = [a.union(b) for a, b in zip(hull, box)]
    # The proved coordinate gap can be shared across all images and all phases.
    margin = ball(Fraction(1, 100))
    counts = []
    for c, e in zip(centers, columns):
        scales = [max(map(abs, column)) for column in e]
        stack = [(full, 0)]
        leaves = depth = 0
        while stack:
            local, level = stack.pop()
            center = enclosed_box(c, e, [(mid, Fraction(0)) for mid, _ in local])
            tangents = [[ball(entry) for entry in edge] for edge in e]
            offsets = [arb(0, ball(radius)) for _, radius in local]
            for _ in range(6):
                enclosure = [
                    center[k] + sum(t[k] * w for t, w in zip(tangents, offsets)) for k in range(3)
                ]
                derivative = step(enclosure, params, tangents)
                if derivative is None:
                    break
                tangents = derivative[1]
                center = step(center, params)[0]
                image = [
                    center[k] + sum(t[k] * w for t, w in zip(tangents, offsets)) for k in range(3)
                ]
                if not any(a > b + margin or a < b - margin for a, b in zip(image, hull)):
                    break
            else:
                depth = max(depth, level)
                leaves += 1
                continue
            if level >= 36:
                raise ArithmeticError("Undecided cell at the subdivision limit")
            j = max(range(3), key=lambda i: scales[i] * local[i][1])
            mid, radius = local[j]
            for sign in (-1, 1):
                child = local.copy()
                child[j] = (mid + sign * radius / 2, radius / 2)
                stack.append((child, level + 1))
        counts.append({"leaves": leaves, "depth": depth})
    return counts


def run():
    ctx.prec = 128
    data = json.loads(INPUT.read_bytes())
    input_hash = sha256(INPUT.read_bytes()).hexdigest()
    for mode, expected_cells in (("primary", 16467), ("ad", 9449)):
        report = INPUT.parent.parent / f"raar-covering-audit-{mode}01/report.json"
        audit = json.loads(report.read_bytes())
        assert audit["complete_graph"] and audit["input_sha256"] == input_hash
        assert audit["total_accepted_cells"] == expected_cells
    assert data["edges"] == [[0, 0], [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 0]]
    centers = [[Fraction(x) for x in c] for c in data["centers"]]
    columns = [
        [tuple(Fraction(row[j]) * Fraction(r[j]) for row in basis) for j in range(3)]
        for basis, r in zip(data["bases"], data["radii"])
    ]
    margins = separated_windows(centers, columns)
    phases = phase_separation(centers, columns, data)

    low, high = Fraction(12554228710, 10**10), Fraction(12554228711, 10**10)
    assert low**7 - low**6 < 1 < high**7 - high**6
    assert ball(low).log() / 7 > ball(Fraction(32496066, 10**9))

    # Trace(A^n) counts itineraries with a marked time origin.
    adjacency = [[int([i, j] in data["edges"]) for j in range(7)] for i in range(7)]
    power = [[int(i == j) for j in range(7)] for i in range(7)]
    traces = [7]
    primitives = {}
    for n in range(1, 101):
        power = [
            [sum(power[i][k] * adjacency[k][j] for k in range(7)) for j in range(7)]
            for i in range(7)
        ]
        traces.append(sum(power[i][i] for i in range(7)))
        if n >= 7:
            assert traces[n] == traces[n - 1] + traces[n - 7]
        primitives[n] = traces[n] - sum(primitives[d] for d in range(1, n) if n % d == 0)
        assert primitives[n] % n == 0
    return {
        "source_window_sha256": input_hash,
        "precision_bits": ctx.prec,
        "separated_window_pairs": margins,
        "phase_gap_strictly_above": "1/100 in the max norm",
        "phase_subdivisions_per_window": phases,
        "perron_root_strict_bracket": [str(low), str(high)],
        "entropy_per_original_step_strictly_above": "0.032496066",
        "orbit_count_lower_bounds_by_least_F_period": {
            str(7 * n): primitives[n] // n for n in (7, 8, 9, 14, 20, 30, 50, 70, 100)
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="save the new receipt")
    args = parser.parse_args()
    result = run()
    if args.write:
        RECEIPT.parent.mkdir(parents=True, exist_ok=True)
        RECEIPT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    else:
        assert result == json.loads(RECEIPT.read_bytes()), "Certificate differs"
    print("All 21 window pairs and all six intermediate phases separated.")
    print(result["phase_subdivisions_per_window"])
    print(result["orbit_count_lower_bounds_by_least_F_period"])

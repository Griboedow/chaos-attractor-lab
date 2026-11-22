"""Isolate the unique cavity equilibrium and its three multipliers with Arb."""

import argparse
import json
from pathlib import Path
from flint import arb, ctx, fmpq
from finite_enclosure import dyadic, NAMES


def run(record):
    ctx.prec = 256
    P, a, d, k, g, h = [dyadic(record["parameters"][k]) for k in NAMES]
    assert d > 0 and k + g > 0 and d + k + g < arb.pi()

    def residual(I):
        phase = d + (k + g) * I / (1 + I)
        return I * (1 + a * a - 2 * a * phase.cos()) - P * P

    low = fmpq(0)
    high = fmpq(2)
    assert residual(arb(low)) < 0 and residual(arb(high)) > 0
    for _ in range(140):
        mid = (low + high) / 2
        v = residual(arb(mid))
        if v > 0:
            high = mid
        elif v < 0:
            low = mid
        else:
            raise ArithmeticError("Increase root precision")
    I = arb(low).union(arb(high))
    T = g * I / (1 + I)
    theta = d + (k + g) * I / (1 + I)
    sn, cs = theta.sin_cos()
    u = a * cs
    v = a * sn
    den = (1 - u) ** 2 + v * v
    x = P * (1 - u) / den
    y = P * v / den
    fx = 2 * x / (1 + I) ** 2
    fy = 2 * y / (1 + I) ** 2
    J = [
        [u - y * k * fx, -v - y * k * fy, -y],
        [v + (x - P) * k * fx, u + (x - P) * k * fy, x - P],
        [h * g * fx, h * g * fy, 1 - h],
    ]

    def characteristic(mu):
        m = [[mu - J[i][j] if i == j else -J[i][j] for j in range(3)] for i in range(3)]
        return (
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        )

    brackets = [("0.05780", "0.05783"), ("1.00159", "1.00162"), ("3.91303", "3.91306")]
    for lo, hi in brackets:
        assert characteristic(arb(lo)) * characteristic(arb(hi)) < 0
    assert arb("6.88810") < T and T < arb("6.88811")
    return dict(
        case=record["case"],
        precision_bits=256,
        bisections=140,
        positive_intensity_rational_bracket=[str(low), str(high)],
        thermal_bracket=["6.88810", "6.88811"],
        coordinates_approx=[float(x), float(y), float(T)],
        multiplier_strict_brackets=brackets,
        stable_dimension=1,
        unstable_dimension=2,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--record", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    r = run(json.loads(a.record.read_text()))
    a.out.write_text(json.dumps(r, indent=2) + "\n")
    print(r)

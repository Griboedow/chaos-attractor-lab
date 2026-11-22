"""Enclose finite paths of both archived initial states by rigorous Arb ball arithmetic."""

import argparse
import json
from pathlib import Path
from flint import arb, ctx, fmpq

NAMES = (
    "pump",
    "field_retention",
    "base_detuning",
    "fast_shift",
    "memory_gain",
    "relaxation",
)


def dyadic(value):
    return arb(fmpq(*float(value).as_integer_ratio()))


def enclose(params, initial, which, bits, start, stop):
    ctx.prec = bits
    P, a, d, k, g, h = map(dyadic, params)
    x, y, t = map(dyadic, initial)
    lower, upper = [5, 6.3] if which == "lower" else [7.6, 11]
    lo = arb(str(lower))
    hi = arb(str(upper))
    limit = arb(fmpq(1, 2**64))
    passed = 0
    for n in range(1, stop + 1):
        I = x * x + y * y
        f = I / (1 + I)
        sn, cs = (d + k * f + t).sin_cos()
        x, y, t = (
            P + a * (cs * x - sn * y),
            a * (sn * x + cs * y),
            (1 - h) * t + h * g * f,
        )
        if not x.is_finite() or not y.is_finite() or not t.is_finite():
            raise ValueError(f"Lost bound at iterate {n}")
        if n >= start:
            if not (lo < t and t < hi):
                raise ValueError(f"Barrier failure at iterate {n}")
            if not all(v.rad() < limit for v in (x, y, t)):
                raise ValueError(f"Enclosure radius reached 2^-64 at iterate {n}")
            passed += 1
    return dict(
        branch=which,
        initial_binary64=initial,
        precision_bits=bits,
        first_iteration=start,
        last_iteration=stop,
        consecutive_steps_certified=passed,
        strict_thermal_bounds=[lower, upper],
        strict_coordinate_radius_limit="2^-64",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--record", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--branch", choices=("lower", "upper"), required=True)
    p.add_argument("--bits", type=int, default=24576)
    p.add_argument("--first", type=int, default=2048)
    p.add_argument("--last", type=int, default=8192)
    a = p.parse_args()
    r = json.loads(a.record.read_text())
    p = [r["parameters"][k] for k in NAMES]
    initial = (
        r["initial"] if a.branch == "lower" else r["other_start_controls"][1]["initial"]
    )
    result = enclose(p, initial, a.branch, a.bits, a.first, a.last)
    a.out.write_text(json.dumps(result, indent=2) + "\n")
    print(result)

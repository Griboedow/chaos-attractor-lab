"""Verify three periodic orbits of cavity Case 56729 and two nearby launch states."""
import argparse
import json
from pathlib import Path

from flint import arb, arb_mat, ctx, fmpq

ROOT = Path(__file__).resolve().parents[3]
RECORD = ROOT / "results/wow_search_20261122/orbit/case-56729.json"
RECEIPT = ROOT / "results/low_hanging_20261122/optics/certificate.json"
IDENTITY = arb_mat([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
SEEDS = {
    "sink": (10, [".9567987", "-.5370079", "-12.9456707"], (".26269", ".26270"), (".8745", ".8749")),
    "gap": (2, [".31330529", ".54460938", "-9.2223509"], ("1.13826", "1.13828"), (".7982", ".7985")),
    "lobes": (2, ["-.09291792", "-.40749198", "-8.52866088"], (".47526", ".47528"), ("1.235", "1.237")),
}
H_VALUES = [["60.8893689", "-27.0516560", "-139.3242328"],
            ["-27.0516560", "20.0440664", "33.8971752"],
            ["-139.3242328", "33.8971752", "450.2685933"]]


def parameters():
    data = json.loads(RECORD.read_text())["parameters"]
    names = ("pump", "field_retention", "base_detuning", "fast_shift", "memory_gain", "relaxation")
    return [arb(fmpq(*data[name].as_integer_ratio())) for name in names]


def ball(center, radius):
    return center + arb(0, radius)


def step(state, p, tangent=False):
    x, y, T = state
    P, a, d, k, g, h = p
    intensity = x*x + y*y
    inverse = 1 / (1 + intensity)
    f = 1 - inverse
    sin, cos = (d + k*f + T).sin_cos()
    u, v = a*(x*cos-y*sin), a*(x*sin+y*cos)
    out = [P+u, v, (1-h)*T+h*g*f]
    if not tangent:
        return out
    fx, fy = 2*x*inverse**2, 2*y*inverse**2
    jac = arb_mat([[a*cos-v*k*fx, -a*sin-v*k*fy, -v],
                   [a*sin+u*k*fx, a*cos+u*k*fy, u],
                   [h*g*fx, h*g*fy, 1-h]])
    return out, jac


def lap(x, n, p):
    J = IDENTITY
    phases = []
    for _ in range(n):
        phases.append(x)
        x, A = centered_step(x, p)
        J = A*J
    return x, J, phases


def centered_step(x, p):
    midpoint = [v.mid() for v in x]
    _, J = step(x, p, tangent=True)
    image = col(step(midpoint, p)) + J*(col(x)-col(midpoint))
    return image.entries(), J


def col(v):
    return arb_mat([[u] for u in v])


def norm(m):
    return max(sum(abs(t) for t in row).upper() for row in m.tolist())


def center(seed, period, p):
    x = [arb(t) for t in seed]
    for _ in range(8):
        end, J, _ = lap(x, period, p)
        correction = (J-IDENTITY).solve(col(end)-col(x))
        x = [(u-correction[i, 0]).mid() for i, u in enumerate(x)]
    return x


def positive_definite(m):
    # Sylvester's criterion for each symmetric matrix contained in the interval matrix.
    return all(arb_mat([row[:n] for row in m.tolist()[:n]]).det() > 0 for n in (1, 2, 3))


def certify(name, x, p):
    n, _, real, mod = SEEDS[name]
    radius = arb("1e-60")
    end, J0, _ = lap(x, n, p)
    A = (J0-IDENTITY).mid().inv().mid()
    box = [ball(t, radius) for t in x]
    _, J, phases = lap(box, n, p)
    q = norm(IDENTITY-A*(J-IDENTITY))
    drift = norm(A*(col(end)-col(x)))
    assert drift + q*radius < radius, (name, float(drift), float(q), float(radius))
    assert all(any(not a.overlaps(b) for a, b in zip(u, v)) for i, u in enumerate(phases) for v in phases[:i]), name
    detuning_ranges = {"sink": ("-12.97","-12.88"), "gap": ("-9.23","-9.17"), "lobes": ("-8.54","-8.46")}
    assert all(labels(t[2], *detuning_ranges[name]) for t in phases)

    # A sign bracket and negative cubic discriminant isolate the real root.
    d, c, b, _ = list(J.charpoly())
    disc = 18*b*c*d - 4*b**3*d + b**2*c**2 - 4*c**3 - 27*d**2
    assert disc < 0
    lo, hi = map(arb, real)
    assert J.charpoly()(lo)*J.charpoly()(hi) < 0
    det = (p[1]**2 * (1-p[5]))**n
    modulus_sq = det / lo.union(hi)
    assert arb(mod[0])**2 < modulus_sq < arb(mod[1])**2, (name, modulus_sq.sqrt())
    return box, dict(period=n, center_approx=[float(t) for t in x], detuning_range=detuning_ranges[name], real_multiplier=real,
                     complex_modulus=mod, newton_row_norm_upper=float(q.upper()) * 1.000001,
                     newton_shift_over_radius_upper=float((drift/radius).upper()) * 1.000001)


def capture_metric(sink, p):
    # The quadratic form was selected by a numerical discrete Lyapunov solve.
    H = arb_mat([[arb(t) for t in row] for row in H_VALUES])
    assert positive_definite(H-IDENTITY)
    r = arb("1e-8")
    _, J, _ = lap([ball(t, 2*r) for t in sink], 10, p)
    assert positive_definite(H-J.transpose()*H*J-IDENTITY*arb(".3"))
    return H, r


def labels(x, lo, hi):
    return arb(lo) < x < arb(hi)


def real_string(t):
    m, e = t.mid().man_exp()
    return [str(m), int(e)]


def restore(dyads):
    return [arb((int(m), e)) for m, e in dyads]


def capture_time(state, sink, H, radius):
    p = parameters()
    for i in range(1, 6001):
        state = step(state, p)
        diff = col([v - c for v, c in zip(state, sink)])
        if norm(diff) < radius and (diff.transpose()*H*diff)[0, 0] < radius**2:
            return i
    raise ArithmeticError("Could not verify entry into the ellipsoid")


def propagate(launches, sink_box, H, capture_radius):
    ctx.prec = 16384
    p = parameters()
    sink = [arb(t) for t in sink_box]
    witnesses = {}
    for sign, dyads in launches.items():
        x = restore(dyads)
        if sign == "minus":
            witnesses[sign] = {"exact_start_as_mantissa_and_binary_exponent": dyads,
                               "sink_ellipsoid_entry_step": capture_time(x, sink, H, capture_radius)}
            continue
        parity = None
        for i in range(1, 4248):
            x = step(x, p)
            if i >= 2200:
                assert labels(x[2], "-8.71", "-8.45"), (i, x[2])
                side = 0 if labels(x[0], "-.46", ".33") else 1 if labels(x[0], "1.88", "2.82") else None
                assert side is not None and (parity is None or side == (parity + i) % 2), (i, x[0])
                if parity is None:
                    parity = (side-i) % 2
        witnesses[sign] = {"exact_start_as_mantissa_and_binary_exponent": dyads,
                           "alternate_lobe_steps_inclusive": [2200,4247], "even_step_lobe": parity}
    return witnesses


def run():
    ctx.prec = 512
    p = parameters()
    centers = {name: center(seed, n, p) for name,(n,seed,*_) in SEEDS.items()}
    roots = {}
    for name,x in centers.items():
        box, roots[name] = certify(name, x, p)
        if name == "sink":
            sink_box = box
    # Keep each initial coordinate an exact dyadic, independent of replay precision.
    v = [arb("-.63427497"), arb("-.59923216"), arb(".48848345")]
    gap = centers["gap"]
    launches = {label:[real_string(t+sign*arb("1e-28")*u) for t,u in zip(gap,v)]
                for label,sign in (("minus",-1),("plus",1))}
    for d in launches.values():
        assert sum((abs(t-c)+arb("1e-60"))**2 for t,c in zip(restore(d),gap)) < arb("1.01e-28")**2
    H, radius = capture_metric(centers["sink"], p)
    witnesses = propagate(launches, sink_box, H, radius)
    starts = json.loads(RECORD.read_text())["other_start_controls"]
    regular = {}
    for index, row in enumerate(starts):
        if row["spectrum"][0] < 0:
            state = [arb(fmpq(*v.as_integer_ratio())) for v in row["initial"]]
            regular[str(index)] = capture_time(state, sink_box, H, radius)
    assert len(regular) == 12
    return {"case":56729,"newton_precision_bits":512,"launch_precision_bits":16384,
            "zero_box_radius":"1e-60","ellipsoid_radius":"1e-8","cycles":roots,
            "witnesses":witnesses,"regular_start_sink_entry_steps":regular}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true",help="Write the verified receipt")
    args = parser.parse_args()
    result = run()
    serial = json.dumps(result, indent=2) + "\n"
    if args.write:
        RECEIPT.write_text(serial,encoding="utf-8",newline="\n")
    else:
        assert json.loads(serial) == json.loads(RECEIPT.read_text())
    print({key: {k:v for k,v in row.items() if k!="exact_start_as_mantissa_and_binary_exponent"}
           for key,row in result["witnesses"].items()})

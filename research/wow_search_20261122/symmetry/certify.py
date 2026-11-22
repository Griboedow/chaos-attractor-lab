"""Verify cyclic trapping annuli in 100-bit Arb from binary64 parameter hex strings."""

import argparse
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx

PACKET = Path(__file__).resolve().parents[3] / "results/wow_search_20261122/symmetry"
ctx.prec = 100


def number(x):
    q = Fraction(x)
    return arb(q.numerator) / q.denominator


from flint import acb

SEEDS = (
    (0.247, 0.193),
    (-0.31, 0.11),
    (0.433, -0.173),
    (0.64, 0.07),
    (-0.216, 0.36),
    (-0.72, -0.08),
    (0.93, 0.04),
)


def certify_starts(p, edges):
    with ctx.workprec(800):
        return entrance_enclosures(p, edges)


def entrance_enclosures(p, edges):
    n = int(p[0])
    l, a, b, g, w, d, e, h = map(number, p[1:])
    hits = []
    for x, y in SEEDS:
        z = acb(number(x), number(y))
        hit = None
        for k in range(160):
            radius = abs(z)
            if radius > 4:
                break
            for j, (lo, hi) in enumerate(edges):
                if radius > lo and radius < hi:
                    hit = {"step": k, "annulus": j}
                    break
            if hit is not None:
                break
            u = z.real**2 + z.imag**2
            z = (
                (acb(l, w) + acb(a, d) * u + b * (z**n).real) * z
                + g * z.conjugate() ** (n - 1)
                + acb(e, h) * z.conjugate() ** (2 * n - 1)
            )
        hits.append({"start": [x, y], "entrance": hit})
    return hits


def forced_graph(pairs):
    ordered = sorted(range(len(pairs)), key=lambda k: Fraction(pairs[k][0]))
    rank = [ordered.index(k) + 1 for k in range(len(pairs))]
    perm = [rank[(i + 1) % len(pairs)] for i in ordered]
    graph = [
        [int(min(left, right) <= k < max(left, right)) for k in range(1, len(pairs))]
        for left, right in zip(perm[:-1], perm[1:])
    ]
    return {
        "radial_ranks_in_time_order": rank,
        "permutation_by_radius": perm,
        "gap_adjacency": graph,
    }


def check(case, pairs, p):
    n = int(p[0])
    l, a, b, g, w, d, e, h = map(number, p[1:])
    edges = [tuple(number(v) for v in pair) for pair in pairs]
    for k, (lo, hi) in enumerate(edges):
        assert lo > 0 and hi > lo
        for v, t in edges[:k]:
            assert (t < lo) or (v > hi)
    # The normalized polar output W is exp(-i theta) F(r exp(i theta)).
    all_rows = []
    for k, (lo, hi) in enumerate(edges):
        nextlo, nexthi = edges[(k + 1) % len(edges)]
        leaves = 0
        maxdepth = 0
        rad_margin = 99.0
        pos_margin = 99.0
        stack = [(lo, hi, arb(0), arb(1), 0)]
        while stack:
            r0, r1, t0, t1, depth = stack.pop()
            r = r0.union(r1)
            s, t = (arb(2) * arb.pi() * t0.union(t1)).sin_cos()
            c2 = t * t - s * s
            s2 = 2 * s * t
            rn = r**n
            rp = rn / r
            r2p = rp * rn
            re = r * (l + a * r * r + b * rn * t) + g * rp * t + r2p * (e * c2 + h * s2)
            im = r * (w + d * r * r) - g * rp * s + r2p * (h * c2 - e * s2)
            norm2 = re**2 + im**2
            # A fixed lower margin for squared radius and positive real polar part.
            margin = number("0.000001")
            lower = norm2 - nextlo**2
            upper = nexthi**2 - norm2
            if lower > margin and upper > margin and re > margin:
                leaves += 1
                rad_margin = min(rad_margin, float(lower.lower()), float(upper.lower()))
                pos_margin = min(pos_margin, float(re.lower()))
                maxdepth = max(maxdepth, depth)
                continue
            if depth > 40:
                raise ValueError(
                    (case, k, depth, str(norm2), str(lower), str(upper), str(re))
                )
            # Binary halves cover the parent, including their common edge.
            # Polar variations scale with powers of r, so radius splits predominate.
            angular_scale = 6.3 * (float(r1.upper()) ** (n - 1)) * (abs(float(g)) + 1)
            if float((r1 - r0).upper()) * 4 >= float((t1 - t0).upper()) * angular_scale:
                mid = (r0 + r1) / 2
                stack.extend(
                    [(r0, mid, t0, t1, depth + 1), (mid, r1, t0, t1, depth + 1)]
                )
            else:
                mid = (t0 + t1) / 2
                stack.extend(
                    [(r0, r1, t0, mid, depth + 1), (r0, r1, mid, t1, depth + 1)]
                )
        all_rows.append(
            dict(
                annulus=k,
                leaves=leaves,
                max_depth=maxdepth,
                approximate_smallest_squared_radius_margin=rad_margin,
                approximate_smallest_real_part=pos_margin,
            )
        )
    return dict(
        binary64_hex=[float(x).hex() for x in p],
        annuli=pairs,
        radial_order=forced_graph(pairs),
        entrance_tests=certify_starts(p, edges),
        inclusive_domain="radius x (n theta / (2 pi) mod 1)",
        bits=ctx.prec,
        strict_proved_squared_radius_and_real_part_margin="0.000001",
        partitions=all_rows,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        type=int,
        choices=[5529, 13745, 205188, 202080],
        nargs="+",
        default=[5529, 13745, 205188, 202080],
    )
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        data = json.loads((PACKET / "annuli.json").read_text())
    else:
        from icons import CACHE, parameter_for

        data = {
            str(i): dict(
                binary64_hex=[x.hex() for x in parameter_for(i)],
                annuli=[
                    [str(a), str(b)]
                    for a, b in json.loads(
                        (CACHE / f"proposed-bands-{i}.json").read_text()
                    )
                ],
            )
            for i in args.case
        }
    out = {}
    for i in args.case:
        item = data[str(i)]
        out[str(i)] = check(
            i, item["annuli"], [float.fromhex(x) for x in item["binary64_hex"]]
        )
        print(i, out[str(i)]["partitions"], flush=True)
        if args.verify:
            assert out[str(i)] == item
    if not args.verify:
        (PACKET / "annuli.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()

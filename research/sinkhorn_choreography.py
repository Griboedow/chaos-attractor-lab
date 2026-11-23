"""Verify attracting three-cluster cycles of fixed-rate Sinkhorn scaling in all dimensions."""

import argparse
import itertools
import json
import math
from pathlib import Path

from flint import arb, arb_mat, ctx, fmpq
from sinkhorn_route import Scaling, as_column, certify, eye, inside_disk

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/sinkhorn-choreography"
SEED = [9.13574344, -3.11542766, -12.25117110, -9.13574344]
DEFAULT_OMEGA = fmpq(7, 4)


class ClusterScaling(Scaling):
    """Log-ratio map on three groups, with the original site multiplicities."""

    def __init__(self, proportions, diagonal=99, omega=DEFAULT_OMEGA):
        self.proportions = proportions
        self.diagonal = arb(diagonal)
        kernel = [
            [v + (self.diagonal if i == j else 0) for j, v in enumerate(proportions)]
            for i in range(3)
        ]
        super().__init__(kernel, arb(omega))
        # Both half steps use the same population weights. The full matrix is
        # symmetric, although the reduced three-row table of group sums is not.
        self.kt = self.k
        self.logb = self.loga

    def own_fraction(self, anchored):
        v = [arb(1)] + [a.exp() for a in anchored]
        mean = sum(p * x for p, x in zip(self.proportions, v))
        return [self.diagonal * x / (mean + self.diagonal * x) for x in v]

    def transverse(self, phases):
        products = [eye(2) for _ in range(3)]
        w, r = self.omega, 1 - self.omega
        for phase in phases:
            new, _ = self.step(phase)
            a = self.own_fraction(phase.entries()[2:])
            b = self.own_fraction(new.entries()[:2])
            for i in range(3):
                block = arb_mat(
                    [[r, -w * a[i]], [-w * b[i] * r, r + w * w * a[i] * b[i]]]
                )
                products[i] = block * products[i]
        return products


def centered_seed(proportions):
    # A numerical continuation through positive weights supplies only a proposal.
    # The final root and all spectral inequalities are independently enclosed.
    state = as_column([arb(x) for x in SEED])
    for k in range(1, 7):
        p = [(1 - arb(k) / 6) / 3 + arb(k) / 6 * x for x in proportions]
        model = ClusterScaling(p)
        for _ in range(4):
            end, jac, _ = model.orbit(state, 3)
            state = (state - (jac - eye(4)).solve(end - state)).mid()
    return [float(x) for x in state.entries()]


def multiplier_upper(matrix):
    if not inside_disk(matrix, arb(1)):
        raise ArithmeticError("A cluster splitting mode was unstable")
    lo, hi = 0, 10000
    while hi - lo > 1:
        mid = (hi + lo) // 2
        if inside_disk(matrix, arb(mid) / 10000):
            hi = mid
        else:
            lo = mid
    return hi / 10000


def verify_profile(sizes):
    n = sum(sizes)
    p = [arb(fmpq(x, n)) for x in sizes]
    model = ClusterScaling(p)
    cycle, phases = certify(model, centered_seed(p), 3)
    for z, leader in zip(phases, (1, 2, 0)):
        x = [arb(0), z[0, 0], z[1, 0]]
        masses = row_masses(model, z)
        if not all(
            i == leader or (x[leader] > x[i] and masses[leader] > masses[i])
            for i in range(3)
        ):
            raise ArithmeticError(
                f"The heaviest individual row did not rotate at {sizes}"
            )
    radii = [multiplier_upper(t) for t in model.transverse(phases)]
    return {
        "sizes": sizes,
        "seed": cycle["seed_approx"],
        "collective_radius_upper": cycle["multiplier_modulus_upper"],
        "splitting_radii_upper": radii,
    }


def row_masses(model, z):
    d = model.n - 1
    u = [arb(1)] + [v.exp() for v in z.entries()[:d]]
    v = [arb(1)] + [v.exp() for v in z.entries()[d:]]
    return [
        ui * sum(kij * vj for kij, vj in zip(row, v)) for ui, row in zip(u, model.k)
    ]


def asymmetric_control():
    n = 6
    # Subtract the Laplacian of a path with pairwise distinct edge weights.
    # All marginals remain exactly one. Every simultaneous permutation
    # symmetry must fix each weighted edge, hence fixes the entire path.
    eps = fmpq(1, 300000)
    kernel = [
        [fmpq(99, 100) * int(i == j) + fmpq(1, 100 * n) for j in range(n)]
        for i in range(n)
    ]
    for i in range(n - 1):
        change = (i + 1) * eps
        kernel[i][i] -= change
        kernel[i + 1][i + 1] -= change
        kernel[i][i + 1] += change
        kernel[i + 1][i] += change
    model = Scaling([[arb(k) for k in row] for row in kernel], arb(7) / 4)
    maximum_radius = 0
    minimum_gap = arb(1)
    seeds = {}
    count_cycles = 0
    for tail in itertools.product(range(3), repeat=n - 1):
        labels = (0,) + tail
        if len(set(labels)) != 3:
            continue
        sizes = tuple(labels.count(k) for k in range(3))
        if sizes not in seeds:
            a, b, c, d = centered_seed([arb(s) / n for s in sizes])
            seeds[sizes] = ([0, a, b], [0, c, d])
        x, y = seeds[sizes]
        seed = [x[k] for k in tail] + [y[k] for k in tail]
        state = as_column([arb(xi) for xi in seed])
        for _ in range(4):
            end, jac, _ = model.orbit(state, 3)
            state = (state - (jac - eye(2 * n - 2)).solve(end - state)).mid()
        cycle, phases = certify(model, [float(q) for q in state.entries()], 3)
        maximum_radius = max(maximum_radius, cycle["multiplier_modulus_upper"])
        for z, leader in zip(phases, (1, 2, 0)):
            masses = row_masses(model, z)
            top_lo = max(q.lower() for q in masses)
            top_hi = max(q.upper() for q in masses)
            for label, q in zip(labels, masses):
                gap = (
                    q.lower() / top_hi - arb(1) / 2
                    if label == leader
                    else arb(1) / 2 - q.upper() / top_lo
                )
                if not gap > 0:
                    raise ArithmeticError(f"Ambiguous heavy-row itinerary at {labels}")
                minimum_gap = min(minimum_gap, gap.lower())
        count_cycles += 1
    if count_cycles != 180:
        raise ArithmeticError("Incorrect control census")
    return {
        "n": n,
        "edge_perturbation": "edge i--i+1: add i/300000, i=1..5; subtract weighted degrees from diagonal",
        "labeled_cycles": count_cycles,
        "max_full_return_radius_upper": maximum_radius,
        "half_maximum_row_threshold_minimum_gap_lower": math.floor(
            1e4 * float(minimum_gap)
        )
        / 1e4,
    }


def canonical(sizes):
    return min(sizes, sizes[1:] + sizes[:1], sizes[2:] + sizes[:2])


def count(n):
    classes = []
    total = 0
    for a in range(1, n - 1):
        for b in range(1, n - a):
            sizes = (a, b, n - a - b)
            if sizes != canonical(sizes):
                continue
            row = verify_profile(sizes)
            # Order of the largest scaling group identifies the cyclic partition.
            labeled = math.factorial(n) // math.prod(math.factorial(k) for k in sizes)
            if len(set(sizes)) == 1:
                labeled //= 3
            row["labeled_cycles"] = labeled
            total += labeled
            classes.append(row)
    expected = (3**n - 3 * 2**n + 3) // 3
    if total != expected:
        raise ArithmeticError(
            "Certified cyclic partitions do not exhaust the declared count"
        )
    return {
        "n": n,
        "labeled_cycles": total,
        "profiles_modulo_relabeling": len(classes),
        "certificates": classes,
    }


def old_control():
    # Equal proportions reduce a .99 example to diagonal 298, other entries 1.
    # The earlier fragile triplet instead has diagonal 8 and relaxation 9/5.
    model = ClusterScaling([arb(1) / 3] * 3, diagonal=fmpq(7, 3), omega=fmpq(9, 5))
    old = json.loads((ROOT / "results/sinkhorn-route/symmetric.json").read_text())
    collective, phases = certify(model, old["cycle"]["seed_approx"], 3)
    for product in model.transverse(phases):
        poly = product.charpoly()
        if not (poly(arb("1.5710")) < 0 and poly(arb("1.5711")) > 0):
            raise ArithmeticError(
                "The old triplet's transverse expanding root is unresolved"
            )
    return {
        "effective_diagonal": 8,
        "omega": "9/5",
        "collective_upper": collective["multiplier_modulus_upper"],
        "splitting_multiplier_interval": [1.5710, 1.5711],
    }


def large_cluster_controls():
    cases = [
        (272, [16.47431, 1.32510, -21.96892, -12.52924], "0.9999248", "0.9999249"),
        (273, [16.48025, 1.32766, -21.97683, -12.53259], "1.0000435", "1.0000436"),
    ]
    records = []
    for n, seed, lo, hi in cases:
        model = ClusterScaling([arb(1) / n, arb(1) / n, arb(n - 2) / n])
        cycle, phases = certify(model, seed, 3)
        poly = model.transverse(phases)[2].charpoly()
        # The 2x2 block has one small real multiplier and a second near +1.
        if not (
            poly(arb(0)) > 0
            and poly(arb("0.2")) < 0
            and poly(arb(lo)) < 0
            and poly(arb(hi)) > 0
        ):
            raise ArithmeticError(
                "Could not separate the large group's two multipliers"
            )
        records.append(
            {
                "n": n,
                "sizes": [1, 1, n - 2],
                "seed": cycle["seed_approx"],
                "collective_radius_upper": cycle["multiplier_modulus_upper"],
                "large_group_leading_multiplier_interval": [lo, hi],
                "multiplicity": n - 3,
                "other_multiplier_interval": [0, 0.2],
            }
        )
    return records


def main(out, upto):
    ctx.prec = 416
    results = []
    for n in range(3, upto + 1):
        row = count(n)
        print(n, row["labeled_cycles"], row["profiles_modulo_relabeling"], flush=True)
        results.append(row)
    out.mkdir(exist_ok=True, parents=True)
    packet = {
        "omega": "7/4",
        "kernel": "0.99 I + 0.01 J/N",
        "marginals": "all row and column sums one",
        "row_roles": [1, 2, 0],
        "old_control": old_control(),
        "asymmetric_control": asymmetric_control(),
        "large_cluster_controls": large_cluster_controls(),
        "dimensions": results,
    }
    (out / "cycles.json").write_text(
        json.dumps(packet, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--upto", type=int, default=18)
    args = parser.parse_args()
    main(args.out, args.upto)

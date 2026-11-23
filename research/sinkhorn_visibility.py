"""Test which certified Sinkhorn cycles are reached from log-random starts."""

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter

import numpy as np
from flint import arb, arb_mat, ctx, fmpq
from sinkhorn_choreography import OUT, canonical, multiplier_upper, row_masses
from sinkhorn_route import Scaling, as_column, certify, enclose_root, eye

N = 18
SEEDS = (271828, 314159, 161803)
SIGMAS = (0.5, 1, 2, 3, 4, 6, 8)
TOL = 1e-7
# A 4-role proposal found from one 18-row Gaussian log start. The first
# three groups lead successively. The fourth never has the largest row.
ORIGINAL = np.array([2, 4, 1, 11]) / 18
START = [
    -13.76710420597,
    -11.41593562168,
    -6.72055355096,
    4.16584684125,
    15.25611620110,
    7.66230686400,
]


def step(x, y, p):
    v = np.exp(y)
    x = -0.75 * x - 1.75 * np.log(0.99 * v + 0.01 * (v * p).sum(axis=-1, keepdims=True))
    u = np.exp(x)
    y = -0.75 * y - 1.75 * np.log(0.99 * u + 0.01 * (u * p).sum(axis=-1, keepdims=True))
    center = (x * p).sum(axis=-1, keepdims=True)
    return x - center, y + center


def four_model(p):
    k = [[arb(99) * int(i == j) + pj for j, pj in enumerate(p)] for i in range(4)]
    model = Scaling(k, arb(7) / 4)
    model.kt, model.logb = model.k, model.loga
    return model


def four_seed(sizes):
    target = np.array(sizes) / sum(sizes)
    state = as_column([arb(x) for x in START])
    for fraction in (0.25, 0.5, 0.75, 1):
        p = ORIGINAL * (1 - fraction) + target * fraction
        m = four_model([arb(float(x)) for x in p])
        for _ in range(5):
            end, j, _ = m.orbit(state, 3)
            state = (state - (j - eye(6)).solve(end - state)).mid()
    return [float(x) for x in state.entries()]


def verify_four(sizes):
    p = [arb(fmpq(x, sum(sizes))) for x in sizes]
    m = four_model(p)
    info, phases = certify(m, four_seed(sizes), 3)
    prods = [eye(2) for _ in p]
    w, r = m.omega, 1 - m.omega
    for leader, z in enumerate(phases):
        weights = row_masses(m, z)
        if not all(i == leader or weights[leader] > q for i, q in enumerate(weights)):
            raise ArithmeticError("Different leading groups")
        new, _ = m.step(z)

        def fractions(logs):
            v = [arb(1)] + [t.exp() for t in logs]
            mean = sum(q * t for q, t in zip(p, v))
            return [99 * t / (99 * t + mean) for t in v]

        for i, (a, b) in enumerate(
            zip(fractions(z.entries()[3:]), fractions(new.entries()[:3]))
        ):
            prods[i] = (
                arb_mat([[r, -w * a], [-w * b * r, r + w * w * a * b]]) * prods[i]
            )
    return {
        "sizes": sizes,
        "seed": info["seed_approx"],
        "collective_upper": info["multiplier_modulus_upper"],
        "splitting_upper": [multiplier_upper(mat) for mat in prods],
    }


def catalogue():
    ctx.prec = 416
    existing = json.loads((OUT / "cycles.json").read_text())["dimensions"][-1][
        "certificates"
    ]
    four, failures = [], []
    for a, b, c in itertools.product(range(1, N - 2), repeat=3):
        d = N - a - b - c
        if d <= 0 or (a, b, c) != canonical((a, b, c)):
            continue
        sizes = (a, b, c, d)
        try:
            four.append(verify_four(sizes))
        except (ArithmeticError, ZeroDivisionError, ValueError) as exc:
            failures.append([sizes, str(exc)])
    equal = verify_four((1, 1, 1, 1))
    return existing, four, failures, equal


def templates(rows):
    result = {}
    for row in rows:
        sizes = np.asarray(row["sizes"])
        p = sizes / sizes.sum()
        h = len(p)
        x = np.r_[0, row["seed"][: h - 1]]
        y = np.r_[0, row["seed"][h - 1 :]]
        x, y = x - x @ p, y + x @ p
        # Relax only the matrix's absolute scale to its unique periodic value.
        for _ in range(90):
            x, y = step(x, y, p)
        phases = []
        leaders = []
        for _ in range(3):
            phases.append(np.stack((x, y)))
            v = np.exp(y)
            leaders.append(np.argmax(x + np.log(0.99 * v + 0.01 * (v @ p))))
            x, y = step(x, y, p)
        np.testing.assert_allclose([x, y], phases[0], atol=1e-10, rtol=0)
        if len(set(leaders)) != 3:
            raise ArithmeticError("Ambiguous orbit")
        spectator = [i for i in range(h) if i not in leaders]
        for t in range(3):
            order = leaders[t:] + leaders[:t] + spectator
            points = np.roll(phases, -t, axis=0)[:, :, order]
            key = tuple(sizes[order])
            if h == 3:
                key = (*key, 0)
            result.setdefault(key, []).append(points)
    return result


def classify(snapshot, refs):
    count = len(snapshot)
    x, y = snapshot[:, :, 0], snapshot[:, :, 1]
    v = np.exp(y)
    mass = x + np.log(0.99 * v + 0.01 * v.mean(axis=2, keepdims=True))
    leaders = mass.max(axis=2, keepdims=True) - mass < TOL * 8
    no_overlap = (leaders.sum(axis=1) <= 1).all(axis=1)
    labels = np.where(leaders.any(axis=1), leaders.argmax(axis=1), 3)
    sizes = np.stack([(labels == i).sum(axis=1) for i in range(4)], axis=1)
    outcome = np.full(count, -1)
    outcome[np.max(abs(snapshot), axis=(1, 2, 3)) < TOL] = 0
    for shape in np.unique(sizes[no_overlap], axis=0):
        key = tuple(shape)
        if key not in refs:
            continue
        ix = np.flatnonzero(no_overlap & (sizes == shape).all(axis=1))
        error = np.full(len(ix), np.inf)
        for ref in refs[key]:
            expected = ref[:, :, labels[ix]].transpose(2, 0, 1, 3)
            error = np.minimum(
                error, np.max(abs(snapshot[ix] - expected), axis=(1, 2, 3))
            )
        outcome[ix[error < TOL]] = 4 if shape[3] else 3
    return outcome


def survey(refs):
    result = []
    p = np.full(N, 1 / N)
    for law in ("Gaussian", "uniform"):
        for sigma in SIGMAS:
            rows = []
            for seed in SEEDS:
                rng = np.random.default_rng(seed)
                shape = (2, 1024, N)
                if law == "Gaussian":
                    start = rng.standard_normal(shape)
                else:
                    start = rng.uniform(-math.sqrt(3), math.sqrt(3), shape)
                x, y = start * sigma
                # Divide out only the gauge: every physical matrix is unchanged.
                center = x.mean(axis=1, keepdims=True)
                x, y = x - center, y + center
                earlier = None
                three = []
                for t in range(1, 901):
                    x, y = step(x, y, p)
                    if t in (448, 449, 450, 898, 899, 900):
                        three.append(np.stack([x, y], axis=1))
                    if t in (450, 900):
                        outcome = classify(np.stack(three, axis=1), refs)
                        three = []
                        if t == 450:
                            earlier = outcome
                counts = Counter(outcome.tolist())
                rows.append(
                    {
                        "seed": seed,
                        "fixed": counts[0],
                        "three_groups": counts[3],
                        "four_groups": counts[4],
                        "unresolved": counts[-1],
                        "changed_since_450": int(np.sum(earlier != outcome)),
                    }
                )
            result.append({"law": law, "sigma": sigma, "batches": rows})
            print(
                law,
                sigma,
                [
                    sum(r[k] for r in rows)
                    for k in (
                        "fixed",
                        "three_groups",
                        "four_groups",
                        "unresolved",
                        "changed_since_450",
                    )
                ],
                flush=True,
            )
    return result


def planted_control(cases, refs):
    rng = np.random.default_rng(1741)
    start = []
    expected = []
    labeled_orbits = []
    for row in cases:
        p = row["sizes"]
        key, points = next(iter(templates([row]).items()))
        labels = rng.permutation(np.repeat(np.arange(4), key))
        point = points[0][0][:, labels] + rng.normal(0, 0.03, (2, N))
        labeled_orbits.append(np.roll(points[0][:, :, labels], -1, axis=0))
        start.append(point)
        expected.append(len(p))
    x, y = np.asarray(start).transpose(1, 0, 2)
    last = []
    for t in range(450):
        x, y = step(x, y, 1 / N)
        if t >= 447:
            last.append(np.stack([x, y], axis=1))
    final = np.stack(last, axis=1)
    np.testing.assert_allclose(final, labeled_orbits, atol=TOL, rtol=0)
    np.testing.assert_array_equal(classify(final, refs), expected)
    return {
        "starts": len(cases),
        "log_noise_sd": 0.03,
        "steps": 450,
        "matched": len(cases),
    }


def balanced_spectators(row):
    """H = transpose after row scaling. Prove H(M) is a cyclic relabeling."""
    model = four_model([arb(fmpq(m, N)) for m in row["sizes"]])
    _, phases = certify(model, row["seed"], 3)
    for order in ([1, 2, 0, 3], [2, 0, 1, 3]):
        perm = arb_mat(
            [
                [int(order[i] == j) - int(order[0] == j) for j in range(1, 4)]
                for i in range(1, 4)
            ]
        )

        def equations(z, perm=perm):
            x, y = z.entries()[:3], z.entries()[3:]
            f, a = model.conditional(y, model.k, model.loga)
            r, w = 1 - model.omega, model.omega
            first = as_column(y) - perm * as_column(x)
            second = as_column(
                [r * t + w * s for t, s in zip(x, f)]
            ) - perm * as_column(y)
            blocks = [[-perm, eye(3)], [r * eye(3), w * a - perm]]
            j = arb_mat(
                6,
                6,
                [
                    blocks[i // 3][k // 3][i % 3, k % 3]
                    for i in range(6)
                    for k in range(6)
                ],
            )
            return as_column(first.entries() + second.entries()), j

        residual, _ = equations(phases[0].mid())
        if max(abs(t).upper() for t in residual.entries()) > arb("1e-25"):
            continue
        half_box, _ = enclose_root(equations, row["seed"])
        center = half_box.mid()
        residual, _ = equations(center)
        _, jbox = equations(half_box)
        tight = center - jbox.solve(residual)
        if all(
            a.contains_interior(b) for a, b in zip(phases[0].entries(), tight.entries())
        ):
            return {
                "sizes": row["sizes"],
                "half_step_relabeling": order,
                "spectator_row_sum_exact": 1,
            }
    raise ArithmeticError("Half-step rotation not enclosed")


def write_certificates():
    existing, four, failures, equal = catalogue()
    cert_count = sum(
        math.factorial(N)
        // math.prod(math.factorial(m) for m in r["sizes"])
        // (3 if len(set(r["sizes"][:3])) == 1 else 1)
        for r in four
    )
    expected = (4**N - 4 * 3**N + 6 * 2**N - 4) // 3
    if failures or cert_count != expected:
        raise ArithmeticError("The four-role census is incomplete")
    print(
        "Four-role profiles:",
        len(four),
        "failed",
        len(failures),
        "cycles",
        cert_count,
        flush=True,
    )
    example = next(row for row in four if tuple(row["sizes"]) == (4, 4, 4, 6))
    packet = {
        "balanced_spectators": balanced_spectators(example),
        "dimension": N,
        "omega": "7/4",
        "kernel": "0.99 I + 0.01 J/N",
        "precision_bits": 416,
        "group_order": "first three are strict per-row mass leaders at phases 0,1,2; fourth never leads",
        "four_group_certificates": four,
        "four_group_certified_count": cert_count,
        "equal_proportions": equal,
    }
    path = OUT / "visibility.certificates.json"
    path.write_text(json.dumps(packet, indent=2) + "\n", encoding="utf-8")
    return existing, four


def main():
    existing, four = write_certificates()
    refs = templates(existing + four)
    packet = {
        "dimension": N,
        "certificate_lf_sha256": hashlib.sha256(
            (OUT / "visibility.certificates.json")
            .read_text(encoding="utf-8")
            .encode("utf-8")
        ).hexdigest(),
        "planted_control": planted_control(existing + four, refs),
        "starts_per_seed_and_amplitude": 1024,
        "steps": 900,
        "second_horizon": 450,
        "all_coordinate_tolerance": TOL,
        "seeds": SEEDS,
        "batches": survey(refs),
    }
    (OUT / "visibility.json").write_text(
        json.dumps(packet, indent=2) + "\n", encoding="utf-8"
    )


def figure():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    certificate = json.loads((OUT / "visibility.certificates.json").read_text())
    sample = json.loads((OUT / "visibility.json").read_text())
    row = next(
        r for r in certificate["four_group_certificates"] if r["sizes"] == [4, 4, 4, 6]
    )
    points = templates([row])[(4, 4, 4, 6)][0]
    p = np.array([4, 4, 4, 6]) / 18
    x, y = points[:, 0], points[:, 1]
    v = np.exp(y)
    mass = np.exp(x) * (0.99 * v + 0.01 * (v @ p)[:, None])
    mass = np.vstack([mass, mass[0]])
    paper, ink, muted = "#faf9f5", "#162c3d", "#4f6271"
    fig = plt.figure(figsize=(7.6, 7.4), facecolor=paper)
    fig.text(
        0.07, 0.945, "The fourth group never leads", size=24, weight="bold", color=ink
    )
    fig.text(
        0.07, 0.902, "One attracting cycle of the 18 × 18 matrix", size=13, color=muted
    )
    ax = fig.add_axes([0.16, 0.575, 0.76, 0.252], facecolor=paper)
    for i, (name, color) in enumerate(
        zip(
            ("A · 4 rows", "B · 4 rows", "C · 4 rows", "D · 6 rows"),
            ("#21776c", "#c5522d", "#7060a5", ink),
        )
    ):
        ax.semilogy(
            range(4),
            mass[:, i],
            label=name,
            color=color,
            marker="o" if i < 3 else "s",
            ls="-" if i < 3 else "--",
            lw=2.3,
            ms=5,
        )
    ax.set_ylabel("Mass per row")
    ax.set_xticks(range(4), ["Step 1", "Step 2", "Step 3", "Step 1"])
    ax.legend(
        ncols=4,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.015),
        frameon=False,
        handlelength=1.3,
        columnspacing=1.05,
        fontsize=9,
    )
    fig.text(
        0.07,
        0.465,
        "Which cycles do random starts reach?",
        size=19,
        weight="bold",
        color=ink,
    )
    fig.text(
        0.07,
        0.428,
        "Fraction reaching any four-group cycle · 3,072 starts per point",
        size=12,
        color=muted,
    )
    low = fig.add_axes([0.16, 0.105, 0.76, 0.255], facecolor=paper)
    for law, color, marker in (
        ("Gaussian", "#21776c", "o"),
        ("uniform", "#c5522d", "s"),
    ):
        rows = [r for r in sample["batches"] if r["law"] == law]
        fractions = [sum(s["four_groups"] for s in r["batches"]) / 3072 for r in rows]
        low.plot(
            [r["sigma"] for r in rows],
            fractions,
            label=law.capitalize(),
            color=color,
            marker=marker,
            lw=2.1,
            ms=6,
            markerfacecolor=paper if law == "uniform" else color,
        )
    low.set_xticks(SIGMAS, [str(x) for x in SIGMAS])
    low.set_yticks([0, 0.5, 1])
    low.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    low.set_xlabel("Standard deviation of initial log weights")
    low.legend(loc="lower right", frameon=False, fontsize=10)
    for axis in (ax, low):
        axis.set_axisbelow(True)
        axis.grid(axis="y", color="#d9deda", lw=0.75)
        axis.spines[["top", "right"]].set_visible(False)
        for name in ("left", "bottom"):
            axis.spines[name].set_color("#a7b0b0")
        axis.tick_params(colors=muted, labelsize=10)
    fig.savefig(OUT / "visibility.png", dpi=145, facecolor=paper)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--certificates-only", action="store_true")
    ap.add_argument("--figure-only", action="store_true")
    args = ap.parse_args()
    if args.certificates_only:
        write_certificates()
    elif args.figure_only:
        figure()
    else:
        main()
        figure()

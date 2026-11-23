"""Certify periodic attractors of the doubly projective Sinkhorn iteration."""

import argparse
import json
from pathlib import Path

from flint import arb, arb_mat, ctx, fmpq

ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "results/sinkhorn-screen"
OUTPUT = ROOT / "results/sinkhorn-route"


def exact(value):
    return arb(fmpq(*float(value).as_integer_ratio()))


def eye(n):
    return arb_mat(n, n, [int(i == j) for i in range(n) for j in range(n)])


def as_column(values):
    return arb_mat(len(values), 1, values)


class Scaling:
    def __init__(self, kernel, omega):
        self.k = kernel
        self.kt = list(zip(*kernel))
        self.n = len(kernel)
        self.omega = omega
        self.loga = [sum(row).log() for row in self.k]
        self.logb = [sum(row).log() for row in self.kt]

    @staticmethod
    def conditional(z, rows, targets):
        ez = [arb(1)] + [u.exp() for u in z]
        weights = [[v * u for v, u in zip(row, ez)] for row in rows]
        sums = [sum(row) for row in weights]
        shifts = [target - s.log() for target, s in zip(targets, sums)]
        conditionals = [[v / s for v in row] for row, s in zip(weights, sums)]
        jac = arb_mat(
            [
                [conditionals[0][j] - row[j] for j in range(1, len(ez))]
                for row in conditionals[1:]
            ]
        )
        return [v - shifts[0] for v in shifts[1:]], jac

    def step(self, z):
        d = self.n - 1
        w = self.omega
        r = 1 - w
        x, y = list(z.entries())[:d], list(z.entries())[d:]
        f, a = self.conditional(y, self.k, self.loga)
        x = [r * u + w * v for u, v in zip(x, f)]
        g, b = self.conditional(x, self.kt, self.logb)
        y = [r * u + w * v for u, v in zip(y, g)]
        ri = eye(d) * r
        blocks = [[ri, a * w], [b * (w * r), ri + (b * a) * w**2]]
        jac = arb_mat(
            2 * d,
            2 * d,
            [
                blocks[i // d][j // d][i % d, j % d]
                for i in range(2 * d)
                for j in range(2 * d)
            ],
        )
        return as_column(x + y), jac

    def orbit(self, z, p):
        jac = eye(z.nrows())
        phases = []
        for _ in range(p):
            phases.append(z)
            z, step_jac = self.step(z)
            jac = step_jac * jac
        return z, jac, phases


def inside_disk(matrix, bound):
    """Real Schur recursion for det(bound*z*I - matrix)."""
    coeff = [a * bound**i for i, a in enumerate(matrix.charpoly().coeffs())]
    while len(coeff) > 1:
        last = coeff[-1]
        first = coeff[0]
        if not (abs(first) < last):
            return False
        # Reciprocal-polynomial Schur transform, scaled to a monic polynomial.
        scale = last * last - first * first
        degree = len(coeff) - 1
        coeff = [
            (last * coeff[i] - first * coeff[degree - i]) / scale
            for i in range(1, degree)
        ] + [arb(1)]
    return True


def enclose_root(evaluate, checkpoint):
    c = as_column([exact(v) for v in checkpoint])
    dimension = c.nrows()
    ident = eye(dimension)
    for _ in range(6):
        residual, j = evaluate(c)
        c = (c - j.solve(residual)).mid()
    # Krawczyk operator on the box X=c+[-r,r]^d.
    radius = arb(2) ** -225
    dx = as_column([arb(0, radius) for _ in range(dimension)])
    box = c + dx
    residual, j = evaluate(c)
    preconditioner = j.inv().mid()
    _, jbox = evaluate(box)
    defect = ident - preconditioner * jbox
    image = -preconditioner * residual + defect * dx
    q = max(
        sum(abs(defect[i, j]) for j in range(dimension)).upper()
        for i in range(dimension)
    )
    if not q < 1 or not all(
        dx[i, 0].contains_interior(image[i, 0]) for i in range(dimension)
    ):
        raise ArithmeticError("root enclosure failed")
    return box, {
        "seed_approx": [float(v) for v in c.entries()],
        "radius_power_of_two": -225,
        "krawczyk_ratio_ceiling": float(
            max(abs(v).upper() for v in image.entries()) / radius
        )
        * 1.001,
    }


def certify(model, checkpoint, period):
    def equations(c):
        end, j, _ = model.orbit(c, period)
        return end - c, j - eye(c.nrows())

    box, info = enclose_root(equations, checkpoint)
    _, jbox, phases = model.orbit(box, period)
    dimension = box.nrows()
    separation = min(
        max(abs(phase[i, 0] - box[i, 0]).lower() for i in range(dimension))
        for phase in phases[1:]
    )
    if not separation > 0:
        raise ArithmeticError("least period unresolved")
    # Search only for a convenient strict upper bound, not a sharp spectral fit.
    if not inside_disk(jbox, arb(1)):
        raise ArithmeticError("attractivity unresolved")
    lo, hi = 0, 10000
    while hi - lo > 1:
        mid = (hi + lo) // 2
        if inside_disk(jbox, arb(mid) / 10000):
            hi = mid
        else:
            lo = mid
    far = max(abs(z).lower() for z in box.entries())
    if not far > 0:
        raise ArithmeticError("could not separate cycle from exact fixed point")
    return dict(
        period=period,
        **info,
        phase_separation_lower=(int(float(separation) * 1000) - 1) / 1000,
        multiplier_modulus_upper=hi / 10000,
    ), phases


def archived(row):
    n, _ = row["shape"]
    case_id = row["case"]
    source = SCREEN / f"extension/n{n}-case{case_id}.json"
    r = json.loads(source.read_text())
    kernel = [[exact(v).exp() for v in line] for line in r["log_kernel"]]
    model = Scaling(kernel, exact(r["omega"]))
    cycle, _ = certify(
        model, r["blocks"][-1]["checkpoint"], row["smallest_tested_recurrence"]
    )
    return dict(shape=n, case=case_id, omega=r["omega"], **cycle)


def expanding_return():
    # Newton search initialized near a close three-step return in the archived
    # final 25,000-point cloud. Only the refined seed is needed for verification.
    source = json.loads((SCREEN / "extension/n4-case64.json").read_text())
    kernel = [[exact(v).exp() for v in line] for line in source["log_kernel"]]
    model = Scaling(kernel, exact(source["omega"]))
    seed = [
        -3.430025365,
        83.14783055,
        123.88685903,
        34.58775169,
        71.82316695,
        -48.62682889,
    ]

    def equations(z):
        end, jac, _ = model.orbit(z, 3)
        return end - z, jac - eye(6)

    box, info = enclose_root(equations, seed)
    _, jac, phases = model.orbit(box, 3)
    separation = min(
        max(abs(phase[i, 0] - box[i, 0]).lower() for i in range(6))
        for phase in phases[1:]
    )
    if not separation > 1:
        raise ArithmeticError("least period unresolved")
    for phase in phases:
        if not inside_disk(model.step(phase)[1], arb(1)):
            raise ArithmeticError("individual phase spectrum unresolved")
    poly = jac.charpoly()
    pairs = [
        (-1670, -1669),
        (-1026, -1025),
        (-935, -934),
        (-261, -260),
        (388, 389),
        (5143, 5144),
    ]
    for lo, hi in pairs:
        if not poly(arb(lo) / 1000) * poly(arb(hi) / 1000) < 0:
            raise ArithmeticError("three-step multiplier enclosure failed")
    return dict(
        shape=4,
        case=64,
        omega=source["omega"],
        period=3,
        **info,
        one_step_spectral_radius_upper=1,
        third_return_real_multipliers=[[lo / 1000, hi / 1000] for lo, hi in pairs],
        phases_approx=[[float(x) for x in z.entries()] for z in phases],
    )


def fold_equations(c):
    # x stores two anchored row logs. P cyclically permutes the three rows
    # and anchors them again. Rotating orbits have y=P*x and T(x,y)=(P^2*x,x).
    x = as_column(c.entries()[:2])
    w = c[2, 0]
    p = arb_mat([[-1, 1], [-1, 0]])
    p2 = p * p
    u = [arb(1)] + [v.exp() for v in x.entries()]
    weights = [[v * (8 if i == j else 1) for j, v in enumerate(u)] for i in range(3)]
    mass = [sum(r) for r in weights]
    prob = [[v / s for v in row[1:]] for row, s in zip(weights, mass)]
    f = as_column([(mass[0] / m).log() for m in mass[1:]])
    a = arb_mat([[prob[0][j] - prob[i][j] for j in range(2)] for i in [1, 2]])

    def cov(row, j, k):
        return row[j] * (int(j == k) - row[k])

    h = [
        [
            [w * (cov(prob[0], j, k) - cov(prob[i], j, k)) for k in range(2)]
            for j in range(2)
        ]
        for i in (1, 2)
    ]
    g = f * w + p2 * x * (1 - w) - p * x
    jac = a * w + p2 * (1 - w) - p
    gw = f - p2 * x
    cof = arb_mat([[jac[1, 1], -jac[1, 0]], [-jac[0, 1], jac[0, 0]]])
    grad_det = [
        sum(cof[i, j] * h[i][j][k] for i in range(2) for j in range(2))
        for k in range(2)
    ]
    grad_det.append(sum(cof[i, j] * (a - p2)[i, j] for i in range(2) for j in range(2)))
    value = as_column(g.entries() + [jac.det()])
    deriv = arb_mat([[jac[i, 0], jac[i, 1], gw[i, 0]] for i in range(2)] + [grad_det])
    # At det=0 these are nonzero left and right null vectors.
    left = [jac[1, 0], -jac[0, 0]]
    right = [jac[0, 1], -jac[0, 0]]
    alpha = sum(left[i] * gw[i, 0] for i in range(2))
    beta = sum(
        left[i] * h[i][j][k] * right[j] * right[k]
        for i in range(2)
        for j in range(2)
        for k in range(2)
    )
    return value, deriv, (jac[0, 0], alpha, beta)


def simple_example():
    model = Scaling(
        [[arb(8 if i == j else 1) for j in range(3)] for i in range(3)], arb(9) / 5
    )
    cycle, phases = certify(model, [2.76, -1.2, -3.96, -2.76], 3)
    # Distinguish the reflected orbit by the largest row sum of diag(u)*K*diag(v).
    itinerary = []
    majority_shares = []
    for phase in phases:
        x = [arb(0), phase[0, 0], phase[1, 0]]
        v = [arb(1), phase[2, 0].exp(), phase[3, 0].exp()]
        rows = [xx.exp() * (sum(v) + 7 * v[i]) for i, xx in enumerate(x)]
        i = next(
            i for i in range(3) if all(i == j or rows[i] > rows[j] for j in range(3))
        )
        itinerary.append(i + 1)
        majority_shares.append(rows[i] / sum(rows))
    if itinerary != [2, 3, 1]:
        raise ArithmeticError("row order unresolved")
    if not all(arb("0.69282") < s < arb("0.69283") for s in majority_shares):
        raise ArithmeticError("largest row share unresolved")
    # An explicit integer start enters a contracting neighborhood of phase 0.
    center = phases[0].mid()
    radius = arb(2) ** -18
    delta = as_column([arb(0, radius)] * 4)
    trap = center + delta
    endpoint, _, _ = model.orbit(center, 6)
    _, deriv, _ = model.orbit(trap, 6)
    image = endpoint - center + deriv * delta
    norm = max(sum(abs(deriv[i, j]) for j in range(4)).upper() for i in range(4))
    if not norm < arb(".858") or not all(
        delta[i, 0].contains_interior(image[i, 0]) for i in range(4)
    ):
        raise ArithmeticError("sixth-return trapping box unresolved")
    seed = as_column([arb(16).log(), -arb(4).log(), -arb(64).log(), -arb(16).log()])
    for _ in range(75):
        seed, _ = model.step(seed)
    if not all(trap[i, 0].contains_interior(seed[i, 0]) for i in range(4)):
        raise ArithmeticError("specified integer start missed the trapping box")
    _, fixed_jac = model.step(as_column([arb(0)] * 4))
    if not inside_disk(fixed_jac, arb(801) / 1000):
        raise ArithmeticError("fixed-point attractivity unresolved")
    fold_box, fold = enclose_root(lambda c: fold_equations(c)[:2], [1.99, -1.66, 1.765])
    pivot, alpha, beta = fold_equations(fold_box)[2]
    if not (abs(pivot) > 0 and alpha * beta < 0):
        raise ArithmeticError("quadratic fold opening toward larger omega unresolved")
    w = fold_box[2, 0]
    lower, upper = arb("1.76514850"), arb("1.76514852")
    if not lower < w or not w < upper:
        raise ArithmeticError("fold parameter out of range")
    fold.update(
        omega_interval=[1.76514850, 1.76514852],
        pivot=float(pivot),
        left_Gw=float(alpha),
        left_Gxx_right_right=float(beta),
    )
    return {
        "kernel_diagonal": 8,
        "kernel_off_diagonal": 1,
        "target_sum": 10,
        "omega": "9/5",
        "attracting_fixed_point": [0] * 4,
        "fixed_modulus_upper": 0.801,
        "cycle": cycle,
        "row_itinerary": itinerary,
        "largest_row_share_interval": [0.69282, 0.69283],
        "starting_rows": [4, 64, 1],
        "starting_columns": [64, 1, 4],
        "steps_to_verified_trap": 75,
        "trap_radius_power_of_two": -18,
        "sixth_return_contraction_upper": 0.858,
        "phase_seeds": [[float(v) for v in phase.entries()] for phase in phases],
        "reflected_by": "swap row labels 2 and 3 and the same column labels",
        "rotating_branch_fold": fold,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=int, help="check one archive ID")
    args = parser.parse_args()
    ctx.prec = 640
    rows = json.loads((SCREEN / "analysis/recurrence.json").read_text())["runs"]
    certified = []
    for row in rows:
        if (row["smallest_tested_recurrence"] or 0) > 1 and args.case in (
            None,
            row["case"],
        ):
            result = archived(row)
            print(
                result["shape"],
                result["case"],
                result["period"],
                result["multiplier_modulus_upper"],
                flush=True,
            )
            certified.append(result)
    if args.case is None:
        OUTPUT.mkdir(exist_ok=True)
        (OUTPUT / "certificates.json").write_text(
            json.dumps(certified, indent=2) + "\n"
        )
        result = simple_example()
        (OUTPUT / "symmetric.json").write_text(json.dumps(result, indent=2) + "\n")
        saddle = expanding_return()
        (OUTPUT / "expanding-return.json").write_text(
            json.dumps(saddle, indent=2) + "\n"
        )
        print(
            "symmetric example",
            result["cycle"]["multiplier_modulus_upper"],
            result["rotating_branch_fold"]["omega_interval"],
        )


if __name__ == "__main__":
    main()

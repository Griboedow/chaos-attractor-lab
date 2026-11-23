"""Enclose the first nonlinear term of the fixed site's resonant return maps."""

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from elliptic_sites import linear_data, sine_series
from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "results/cubic-elliptic-sites"
# Decimal proposals. The verifier isolates each root by strict endpoint signs.
PROPOSALS = [
    (3, "0.133642692140877899299454771098939107566828293588971854244795020125"),
    (3, "0.141705738257650019153077162544862259873915293885493670379171539909"),
    (3, "0.142626306990232678640550185847502819228788044749912188818162968374"),
    (4, "0.133634498338077896718395176143994489106215358996015633858942281885"),
    (4, "0.141734706307067891265590101791627015338941410876881270293304985770"),
    (4, "0.142603078096242743887199751238487165243419363315151892847554549185"),
]


def generators(degree):
    powers = [(n - j, j) for n in range(1, degree + 1) for j in range(n + 1)]
    # Re and Im of (a+ib)(conj(h)+conj(h)^2/2) -18*k*h*|h|^2.
    fields = [
        ({(1, 0): 1, (2, 0): 0.5, (0, 2): -0.5}, {(0, 1): -1, (1, 1): -1}),
        ({(0, 1): 1, (1, 1): 1}, {(1, 0): 1, (2, 0): 0.5, (0, 2): -0.5}),
        ({(3, 0): -18, (1, 2): -18}, {(2, 1): -18, (0, 3): -18}),
    ]
    matrices = []
    for components in fields:
        mat = arb_mat(len(powers), len(powers))
        for i, (p, q) in enumerate(powers):
            for axis, (power, component) in enumerate(zip((p, q), components)):
                if not power:
                    continue
                for (u, v), weight in component.items():
                    key = (p + u - (axis == 0), q + v - (axis == 1))
                    if key in powers:
                        j = powers.index(key)
                        mat[i, j] = mat[i, j] + power * weight
        matrices.append(mat)
    return powers, matrices


def jet(gain, degree, steps, order):
    powers, (U, V, W) = generators(degree)
    m = len(powers)
    zero = arb_mat(m, m)
    identity = arb_mat([[int(i == j) for j in range(m)] for i in range(m)])
    h = 2 * arb.pi() / steps
    R = arb(8)
    size = R * h
    ab = abs(gain) * (
        12 + arb(3) / 2 * size.cosh() + 24 * (3 * size).cosh() + 18 * (6 * size).cosh()
    )
    # On the complex time disk ab bounds |a|+|b|. Bound every row of the
    # lifted generator in the induced infinity norm.
    norm = max(
        (
            max(
                sum(abs(U[i, j]) for j in range(m)), sum(abs(V[i, j]) for j in range(m))
            )
            * ab
            + sum(abs(W[i, j]) for j in range(m)) * abs(gain)
        ).upper()
        for i in range(m)
    )
    tail = (size * norm).exp() * R ** (-order) / (R - 1)
    error = arb(0, tail.upper())
    total = identity
    for cell in range(steps):
        s1, s3, s6 = [sine_series(j, cell * h, h, order) for j in (1, 3, 6)]
        coefficients = []
        for i in range(order):
            a = gain * ((12 if i == 0 else 0) - 3 * s1[i] / 2 + 24 * s3[i])
            b = -18 * gain * s6[i]
            coefficients.append(U * a + V * b + (W * gain if i == 0 else zero))
        terms = [identity]
        for n in range(order):
            conv = sum((coefficients[i] * terms[n - i] for i in range(n + 1)), zero)
            terms.append(conv * (h / (n + 1)))
        local = sum(terms, zero)
        total = (
            arb_mat(
                m,
                m,
                [
                    local[i, j] + error if sum(powers[j]) >= sum(powers[i]) else arb(0)
                    for i in range(m)
                    for j in range(m)
                ],
            )
            * total
        )
    return powers, total


def verify(job):
    turns, proposal, steps, order = job
    ctx.prec = 352
    mid = arb(proposal)
    eps = arb("1e-60")
    signs = []
    for endpoint in (mid - eps, mid + eps):
        m, _ = linear_data(endpoint, 1, 128, 94, False)
        signs.append(m[0, 0] + m[1, 1] + int(turns == 3))
    if not signs[0] * signs[1] < 0:
        raise ArithmeticError("Resonance isolation failed")
    gain = mid + arb(0, eps.upper())
    degree = turns - 1
    powers, J = jet(gain, degree, steps, order)
    repeated = J**turns
    cosine = arb(-1) / 2 if turns == 3 else arb(0)
    sine = (1 - cosine**2).sqrt()
    shear = (J[0, 0] - cosine) / sine
    height = J[1, 0] / sine

    def displacement(x, y):
        X, Y = x + shear * y, height * y
        terms = [X**p * Y**q if p + q == degree else arb(0) for p, q in powers]
        re, im = [
            sum(repeated[i, j] * v for j, v in enumerate(terms)) for i in range(2)
        ]
        return acb(re - shear * im / height, im / height)

    first = displacement(arb(1), arb(0))
    packet = {
        "turns": turns,
        "gain_center": proposal,
        "gain_radius": "1e-60",
        "trace_endpoint_residuals": [t.str(8) for t in signs],
        "frame_second_column": [shear.str(16), height.str(16)],
    }
    if turns == 3:
        if not abs(first) > 1:
            raise ArithmeticError(f"Quadratic coefficient unresolved: {first}")
        packet.update(
            coefficient=[first.real.str(20), first.imag.str(20)],
            coefficient_modulus=abs(first).str(16),
            conclusion="unstable",
        )
    else:
        root = 1 / arb(2).sqrt()
        diagonal = acb(root, -root) * displacement(root, root)
        A, B = (first + diagonal) / 2, (first - diagonal) / 2
        # Liouville gives the same real coefficient without the cancellation
        # between cubic terms in four successive coordinate transformations.
        mon, gram = linear_data(gain, 1, 128, 60)
        v = arb_mat([[1], [0]])
        energy = arb(0)
        for _ in range(4):
            energy += (v.transpose() * gram * v)[0, 0]
            v = mon * v
        real = -18 * gain * energy
        if not A.real.overlaps(real) or not real < 0:
            raise ArithmeticError("The independent divergence coefficient differs")
        A = acb(real, A.imag)
        gap = abs(A) - abs(B)
        if gap > 0 and A.real < 0:
            status = "attracting"
        elif gap < 0:
            status = "unstable"
        else:
            raise ArithmeticError(f"Quartic test unresolved: {A}; {B}; {gap}")
        packet.update(
            cubic_isotropic=[A.real.str(20), A.imag.str(20)],
            cubic_resonant=[B.real.str(20), B.imag.str(20)],
            modulus_gap=gap.str(16),
            conclusion=status,
        )
    print(turns, proposal, packet, flush=True)
    return packet


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=128)
    parser.add_argument("--order", type=int, default=112)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--out", type=Path, default=OUT / "resonances.json")
    args = parser.parse_args()
    with ProcessPoolExecutor(args.workers) as pool:
        rows = list(
            pool.map(verify, [(t, k, args.steps, args.order) for t, k in PROPOSALS])
        )
    args.out.write_text(
        json.dumps(
            {
                "case": 48000001,
                "bits": 352,
                "steps": args.steps,
                "order": args.order,
                "site": 1,
                "resonances": rows,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

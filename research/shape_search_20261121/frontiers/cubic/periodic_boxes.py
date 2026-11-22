"""Verify periodic orbits by interval multiple shooting. The only input is a node proposal."""

import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from math import factorial
from pathlib import Path

import numpy as np
from flint import acb, acb_poly, arb, arb_mat, ctx


def coefficient_norm(polynomial):
    return sum((abs(c).upper() for c in polynomial.coeffs()), arb(0)).upper()


def conjugate(polynomial):
    return acb_poly([c.conjugate() for c in polynomial.coeffs()])


def convolution(left, right, degree):
    return sum((left[k] * right[degree - k] for k in range(degree + 1)), acb(0))


def drive_series(theta, step, degree):
    coeffs = [acb(0) for _ in range(degree + 1)]
    remainder = arb(0)
    for frequency, amplitude in [(1, arb(1) / 4), (3, -4), (6, acb(0, 3))]:
        angle, h = frequency * theta, frequency * step
        values = [angle.sin(), h * angle.cos()]
        for k in range(2, degree + 1):
            values.append(-h * h * values[-2] / (k * (k - 1)))
        for k, value in enumerate(values):
            coeffs[k] += amplitude * value
        remainder += abs(amplitude) * h ** (degree + 1) / factorial(degree + 1)
    return coeffs, remainder.upper()


def local_polynomials(start, columns, b, step, gain):
    # Taylor polynomials in u=(theta-theta_0)/step for a state and two tangent columns.
    z = [start]
    v = [[column] for column in columns]
    square, derivative, conjugate_derivative, residual = [], [], [], []
    linear, antilinear = [], []
    for k in range(len(b) - 1):
        square.append(convolution(z, z, k))
        derivative.append(3 * square[k] - (3 if k == 0 else 0))
        conjugate_derivative.append(derivative[k].conjugate())
        residual.append(convolution(square, z, k) - 3 * z[k] + b[k])
        linear.append(-gain * convolution(derivative, conjugate_derivative, k))
        antilinear.append(
            -6 * gain * sum(z[j].conjugate() * residual[k - j] for j in range(k + 1))
        )
        z.append(
            -step * gain * convolution(conjugate_derivative, residual, k) / (k + 1)
        )
        for col in v:
            rhs = convolution(linear, col, k)
            rhs += sum(antilinear[j] * col[k - j].conjugate() for j in range(k + 1))
            col.append(step * rhs / (k + 1))
    return acb_poly(z), [acb_poly(column) for column in v]


def cell(start, columns, errors, forcing, step, gain):
    # Euclidean error disks for the nominal state, the initial box, and both derivative columns.
    nominal_error, box_error, column_errors = errors
    b, tail = forcing
    z, vs = local_polynomials(start, columns, b, step, gain)
    dp = 3 * z * z - 3
    p = z * z * z - 3 * z + acb_poly(b)
    dp_bar = conjugate(dp)
    defect = coefficient_norm(z.derivative() + step * gain * dp_bar * p)
    defect += step * gain * coefficient_norm(dp) * tail

    cap = arb(1) / 1024
    span = coefficient_norm(z - acb_poly([start])) + cap
    size = abs(start) + span
    initial_dp = abs(3 * start * start - 3)
    dp_change = 6 * abs(start) * span + 3 * span**2
    dp_max = initial_dp + dp_change
    dp_min = max(arb(0), (initial_dp - dp_change).lower())
    p_max = (
        abs(start**3 - 3 * start + b[0])
        + initial_dp * span
        + 3 * abs(start) * span**2
        + span**3
    )
    p_max += step * 121 / 4  # uniform bound for |b'|
    mu = (gain * (-(dp_min**2) + 6 * size * p_max)).upper()
    growth = (step * mu).exp().upper()
    peak = max(arb(1), growth)
    box_peak = (peak * (box_error + defect)).upper()
    if not box_peak < cap:
        raise ArithmeticError("State error left the logarithmic-norm enclosure")
    jacobian_lipschitz = gain * (18 * size * dp_max + 6 * p_max)

    a = -gain * dp * dp_bar
    c = -6 * gain * conjugate(z) * p
    next_columns, next_errors = [], []
    for v, error in zip(vs, column_errors):
        v_size = coefficient_norm(v)
        v_defect = coefficient_norm(v.derivative() - step * (a * v + c * conjugate(v)))
        v_defect += (
            step * (6 * gain * size * tail + jacobian_lipschitz * box_peak) * v_size
        )
        endpoint = v(1)
        center = endpoint.mid()
        next_columns.append(center)
        next_errors.append(
            (growth * error + peak * v_defect + abs(endpoint - center)).upper()
        )
    endpoint = z(1)
    center = endpoint.mid()
    rounding = abs(endpoint - center)
    nominal_error = (growth * nominal_error + peak * defect + rounding).upper()
    box_error = (growth * box_error + peak * defect + rounding).upper()
    return center, next_columns, (nominal_error, box_error, next_errors)


def verify_orbit(job):
    index, nodes, subdivisions, degree, bits = job
    ctx.prec = bits
    gain = arb(17) / 50
    radius = arb(2) ** -40
    count = len(nodes)
    step = 2 * arb.pi() / (count * subdivisions)
    dim = 2 * count
    derivative = arb_mat(dim, dim)
    mismatch = arb_mat(dim, 1)
    monodromy = arb_mat([[1, 0], [0, 1]])
    for i, node in enumerate(nodes):
        start = acb(float(node.real), float(node.imag))
        columns = [acb(1), acb(0, 1)]
        errors = (arb(0), arb(2).sqrt() * radius, [arb(0), arb(0)])
        for j in range(subdivisions):
            theta = (i * subdivisions + j) * step
            forcing = drive_series(theta, step, degree)
            start, columns, errors = cell(start, columns, errors, forcing, step, gain)
        following = (i + 1) % count
        target = nodes[following]
        mismatch[2 * i, 0] = arb(float(target.real)) - start.real + arb(0, errors[0])
        mismatch[2 * i + 1, 0] = (
            arb(float(target.imag)) - start.imag + arb(0, errors[0])
        )
        block = arb_mat(2, 2)
        for j, (column, error) in enumerate(zip(columns, errors[2])):
            for k, value in enumerate([column.real, column.imag]):
                block[k, j] = value + arb(0, error)
                derivative[2 * i + k, 2 * i + j] = -block[k, j]
        derivative[2 * i, 2 * following] += 1
        derivative[2 * i + 1, 2 * following + 1] += 1
        monodromy = block * monodromy

    # The numerical inverse is only a proposal. All subsequent tests use Arb.
    midpoint = np.array(
        [[float(derivative[i, j].mid()) for j in range(dim)] for i in range(dim)]
    )
    inverse = arb_mat(np.linalg.inv(midpoint).tolist())
    identity = arb_mat(np.eye(dim, dtype=int).tolist())
    error_matrix = identity - inverse * derivative
    q = max(
        sum(abs(error_matrix[i, j]).upper() for j in range(dim)).upper()
        for i in range(dim)
    )
    offset = max(abs(v).upper() for v in inverse * mismatch)
    if not offset + q * radius < radius:
        raise ArithmeticError(
            f"Orbit {index}: interval Newton inclusion failed, q={q}, offset={offset}"
        )
    # Liouville's formula and div(f) = -2*kappa*|3z^2-3|^2 give 0 < det <= 1.
    trace = monodromy[0, 0] + monodromy[1, 1]
    if trace > 2:
        kind = "positive saddle"
    elif trace < -2:
        kind = "negative saddle"
    else:
        norm = max(
            sum(abs(monodromy[i, j]) for j in range(2)).upper() for i in range(2)
        )
        if not norm < 1:
            raise ArithmeticError(f"Orbit {index}: stability undecided")
        kind = "sink"
    print(f"{index:02d}  {kind}  q={float(q):.4g}", flush=True)
    return {
        "orbit": index,
        "type": kind,
        "contraction_bound": str(q),
        "center_correction_bound": str(offset),
        "return_trace_enclosure": str(trace),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nodes", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--orbit", type=int, nargs="+")
    parser.add_argument("--subdivisions", type=int, default=16)
    parser.add_argument("--degree", type=int, default=14)
    args = parser.parse_args()
    nodes = np.load(args.nodes)["nodes"]
    chosen = list(range(len(nodes))) if args.orbit is None else args.orbit
    bits = 112
    ctx.prec = bits
    for n, index in enumerate(chosen):
        for other in chosen[:n]:
            # At least one node coordinate separates the boxes at a common phase.
            first, second = nodes[index].view(float), nodes[other].view(float)
            if not any(
                abs(arb(float(a)) - arb(float(b))) > arb(2) ** -39
                for a, b in zip(first, second)
            ):
                raise ArithmeticError("Overlapping orbit proposals")
    jobs = [(i, nodes[i], args.subdivisions, args.degree, bits) for i in chosen]
    if args.workers == 1:
        records = list(map(verify_orbit, jobs))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            records = list(pool.map(verify_orbit, jobs))
    info = {
        "gain": [17, 50],
        "case": 48000001,
        "phase_nodes_per_turn": nodes.shape[1],
        "box_half_width": "2^-40",
        "subdivisions": args.subdivisions,
        "degree": args.degree,
        "arb_bits": bits,
        "node_file_sha256": hashlib.sha256(args.nodes.read_bytes()).hexdigest(),
        "orbits": records,
    }
    args.out.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")

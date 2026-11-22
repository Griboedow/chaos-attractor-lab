"""Enclose the linear return and leading radial decay at a stationary cubic site."""

import argparse
import json
from pathlib import Path

from flint import arb, arb_mat, arb_poly, ctx, fmpq

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "results/cubic-elliptic-sites"


def eye():
    return arb_mat([[1, 0], [0, 1]])


def sine_series(freq, t, h, degree):
    angle, step = freq * t, freq * h
    terms = [angle.sin(), step * angle.cos()]
    for j in range(2, degree):
        terms.append(-terms[-2] * step**2 / (j * (j - 1)))
    return terms


def widened(matrix, radius):
    return arb_mat(2, 2, [v + arb(0, radius.upper()) for v in matrix.entries()])


def local_flow(t, h, kappa, site, degree, tail_f, tail_g, need_gram):
    s1, s3, s6 = [sine_series(j, t, h, degree) for j in (1, 3, 6)]
    coefficients = []
    for j in range(degree):
        re = kappa * ((12 if j == 0 else 0) - 6 * site * (s1[j] / 4 - 4 * s3[j]))
        im = -18 * kappa * site * s6[j]
        coefficients.append(arb_mat([[re, im], [im, -re]]))
    f = [eye()]
    gram = arb_mat(2, 2)
    for j in range(degree):
        next_f = sum((coefficients[i] * f[j - i] for i in range(j + 1)), arb_mat(2, 2))
        if need_gram:
            next_g = sum(
                (f[i].transpose() * f[j - i] for i in range(j + 1)), arb_mat(2, 2)
            )
            gram += next_g * (h / (j + 1))
        f.append(next_f * (h / (j + 1)))
    return widened(sum(f, arb_mat(2, 2)), tail_f), widened(gram, tail_g)


def linear_data(kappa, site, n, degree, need_gram=True):
    h = 2 * arb.pi() / n
    # Complexify the two real variational coordinates. On |u| <= R, ||A(t+h*u)||∞ <= M.
    # F and G'=h F^T F are analytic, with entry bounds E and 2 R h E² respectively
    # for F and the integral G. Cauchy's coefficient estimate bounds their finite tails.
    R = arb(8)
    size = R * h
    M = abs(kappa) * (
        12 + arb(3) / 2 * size.cosh() + 24 * (3 * size).cosh() + 18 * (6 * size).cosh()
    )
    E = (size * M).exp()
    geom = R ** (-degree) / (R - 1)
    tail_f = E * geom
    tail_g = 2 * size * E**2 * geom
    phi, gram = eye(), arb_mat(2, 2)
    for j in range(n):
        step, integral = local_flow(
            j * h, h, kappa, site, degree, tail_f, tail_g, need_gram
        )
        if need_gram:
            gram += phi.transpose() * integral * phi
        phi = step * phi
    return phi, gram


def analyze(kappa, site, n, degree):
    phi, gram = linear_data(kappa, site, n, degree)
    a, b, c, d = phi.entries()
    trace = a + d
    if not (trace > -2 and trace < 2 and abs(trace + 1) > 0 and abs(trace) > 0):
        raise ArithmeticError("Nonresonant ellipticity through degree three unresolved")
    q = (1 - (trace / 2) ** 2).sqrt()
    if b < 0:
        q = -q
    elif not b > 0:
        raise ArithmeticError("Rotation orientation unresolved")
    H = arb_mat([[b, (d - a) / 2], [(d - a) / 2, -c]]) * (1 / q)
    # det(phi) = 1 exactly by the trace-free variational equation. This formula
    # gives H positive with det H=1 and phi H phi^T=H, so S=H^(1/2) is area-normalized.
    product = H * gram
    gamma = 9 * kappa * (product[0, 0] + product[1, 1])
    amplitude = 1 / (2 * gamma).sqrt()
    for enclosure, lo, hi in [
        (trace, "-1.20400", "-1.20399"),
        (gamma, "5120.408", "5120.410"),
        (amplitude, "0.009881", "0.009883"),
    ]:
        if not arb(lo) < enclosure or not enclosure < arb(hi):
            raise ArithmeticError(f"Bound failed: {enclosure}, expected ({lo}, {hi})")
    matrix_strings = lambda A: [[A[i, j].str(24) for j in range(2)] for i in range(2)]
    return {
        "site": site,
        "trace": trace.str(24),
        "gamma": gamma.str(24),
        "limiting_scaled_radius": amplitude.str(24),
        "monodromy": matrix_strings(phi),
        "covariance": matrix_strings(H),
        "integrated_gramian": matrix_strings(gram),
    }


def window(center, bracket):
    # Interpolate trace as an analytic function of gain on x=(gain-center)/h.
    # The contour estimate bounds both the approximation and its derivative
    # on the entire real interval; it does not infer monotonicity from samples.
    center = arb(center)
    h = arb(fmpq(1, 10000))
    degree = 18
    points = [arb(fmpq(2 * i - degree, degree)) for i in range(degree + 1)]
    diffs = []
    for x in points:
        p, _ = linear_data(center + h * x, 1, 96, 36, need_gram=False)
        diffs.append(p[0, 0] + p[1, 1])
    interp = arb_poly([diffs[0]])
    basis = arb_poly([1])
    for order in range(1, degree + 1):
        diffs = [
            (diffs[j + 1] - diffs[j]) / (points[j + order] - points[j])
            for j in range(len(diffs) - 1)
        ]
        basis *= arb_poly([-points[order - 1], 1])
        interp += basis * diffs[0]
    R = arb(360)
    # ||B(t)||_infty <= 55.5 at real times for A=kappa*B. Hence the trace obeys
    # |trace Phi(T)| <= 2*exp(T*(|center|+h R)*55.5) on the complex gain contour.
    M = 2 * (2 * arb.pi() * (abs(center) + h * R) * arb("55.5")).exp()
    factor = R * M / (R - 1) ** (degree + 1)
    tail = factor * 2 ** (degree + 1) / (R - 1)
    deriv_tail = factor * (
        (degree + 1) * 2**degree / (R - 1) + 2 ** (degree + 1) / (R - 1) ** 2
    )
    slope = interp.derivative()(arb(0, 1)) + arb(0, deriv_tail.upper())
    ascending = slope > 0
    if not ascending and not slope < 0:
        raise ArithmeticError("Trace monotonicity unresolved")
    values = (-2, 2) if ascending else (2, -2)
    for (lo, hi), level in zip(bracket, values):
        left = interp((arb(lo) - center) / h) + arb(0, tail.upper()) - level
        right = interp((arb(hi) - center) / h) + arb(0, tail.upper()) - level
        if not (left < 0 and right > 0 if ascending else left > 0 and right < 0):
            raise ArithmeticError(
                f"Boundary {level} unresolved at {lo}, {hi}: {left}, {right}"
            )
    return {
        "center": center.str(10, radius=False),
        "gain_half_width": "0.0001",
        "gain_boundary_brackets": bracket,
        "trace_at_boundaries": values,
        "slope_per_gain": [(slope / h).lower().str(12), (slope / h).upper().str(12)],
        "parameter_degree": degree,
        "time_steps": 96,
        "time_degree": 36,
        "trace_interpolation_error": tail.str(8),
        "trace_interpolation_derivative_error": (deriv_tail / h).str(8),
        "polynomial_ascending_coefficients": [v.str(20) for v in interp.coeffs()],
    }


def draw(packet, path):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D

    paper, ink, muted, accent = "#fbfaf5", "#223844", "#526675", "#137d75"
    fig = plt.figure(figsize=(6.5, 6.7), facecolor=paper)
    fig.text(
        0.065,
        0.93,
        "Three narrow windows of attraction",
        fontsize=20.3,
        weight="bold",
        color=ink,
    )
    fig.text(
        0.067,
        0.871,
        "Both stationary points of the forced cubic,  z = ±1",
        fontsize=13,
        color=muted,
    )
    axes = []
    for i, w in enumerate(packet["windows"]):
        center, span = float(w["center"]), float(w["gain_half_width"])
        ax = fig.add_axes([0.135, 0.622 - 0.234 * i, 0.80, 0.162], facecolor=paper)
        axes.append(ax)
        coeff = [float(arb(v)) for v in w["polynomial_ascending_coefficients"]]
        x = np.linspace(-1, 1, 1201)
        trace = np.polynomial.polynomial.polyval(x, coeff)
        bound = [
            np.mean([float(v) for v in pair]) for pair in w["gain_boundary_brackets"]
        ]
        ax.axvspan(*bound, color="#dcefe6", zorder=0)
        for level in (-2, 2):
            ax.axhline(
                level, color="#acbcb4", linewidth=0.85, linestyle=(0, (3, 3)), zorder=1
            )
        ax.plot(center + span * x, trace, color=ink, linewidth=1.7)
        valid = np.abs(trace) <= 2
        ax.plot((center + span * x)[valid], trace[valid], color=accent, linewidth=3.3)
        for level in (-2, -1, 0, 2):
            # The local trace is strictly monotone across the whole displayed interval.
            tr, xx = (trace, x) if trace[0] < trace[-1] else (trace[::-1], x[::-1])
            gain = center + span * np.interp(level, tr, xx)
            ax.scatter(
                gain,
                level,
                s=32,
                facecolor=paper,
                edgecolor=accent,
                linewidth=1.3,
                zorder=4,
            )
        ax.set(xlim=(center - span, center + span), ylim=(-3.4, 3.4), yticks=[-2, 2])
        xticks = center + span * np.array([-0.8, 0, 0.8])
        ax.set_xticks(xticks, labels=[f"{v:.6f}" for v in xticks])
        ax.tick_params(colors=muted, labelsize=9, width=0.7, length=3)
        ax.spines[["right", "top", "left"]].set_visible(False)
        ax.spines["bottom"].set_color("#90a39c")
        if i == 1:
            dot = 0.1417
            height = np.polynomial.polynomial.polyval((dot - center) / span, coeff)
            ax.scatter(dot, height, s=34, color="#984e83", zorder=5)
            ax.annotate(
                "κ = 0.1417",
                xy=(dot, height),
                xytext=(-44, 19),
                textcoords="offset points",
                fontsize=9.5,
                color="#803f6b",
                ha="right",
                arrowprops={"arrowstyle": "-", "color": "#9d8296", "linewidth": 0.8},
            )
    fig.text(
        0.035,
        0.455,
        "Trace of the one-period linear return",
        rotation=90,
        fontsize=10,
        color=muted,
        va="center",
        ha="center",
    )
    axes[-1].set_xlabel("Gain κ", fontsize=11, color=muted, labelpad=4)
    fig.legend(
        handles=[
            Line2D([], [], color=accent, lw=3, label="Algebraic attraction"),
            Line2D([], [], color=ink, lw=1.7, label="Saddle"),
            Line2D(
                [],
                [],
                marker="o",
                linestyle="none",
                markeredgecolor=accent,
                markerfacecolor=paper,
                label="Resonance left open",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.50, 0.857),
        ncol=3,
        frameon=False,
        fontsize=9.6,
        columnspacing=1.1,
        handlelength=1.3,
        labelcolor=ink,
    )
    fig.text(
        0.067,
        0.025,
        "Each panel spans 0.0002 in gain. Window boundaries are interval-verified.",
        color=muted,
        fontsize=10,
    )
    fig.savefig(path, dpi=165, facecolor=paper)
    plt.close(fig)


def main(out):
    ctx.prec = 208
    gain = arb(fmpq(1417, 10000))
    n, degree = 128, 48
    records = [analyze(gain, site, n, degree) for site in (-1, 1)]
    windows = [
        window(
            "0.133635", [["0.13361815", "0.13361816"], ["0.13365090", "0.13365091"]]
        ),
        window(
            "0.141735", [["0.14167810", "0.14167812"], ["0.14179777", "0.14179780"]]
        ),
        window(
            "0.142603", [["0.14255143", "0.14255145"], ["0.14264817", "0.14264820"]]
        ),
    ]
    out.mkdir(exist_ok=True, parents=True)
    packet = {
        "case": 48000001,
        "gain": "1417/10000",
        "steps": n,
        "degree": degree,
        "arb_bits": ctx.prec,
        "period": "2*pi",
        "certified": records,
        "windows": windows,
    }
    (out / "weak-foci.json").write_text(
        json.dumps(packet, indent=2) + "\n", encoding="utf-8"
    )
    draw(packet, out / "windows.png")
    for r in records:
        print(r["site"], r["trace"], r["gamma"], r["limiting_scaled_radius"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    main(parser.parse_args().out)

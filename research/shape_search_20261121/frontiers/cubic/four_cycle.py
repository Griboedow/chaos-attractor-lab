"""Shoot both period-four satellites of the last 1:4 resonance in the forced cubic."""

import argparse
import json
from pathlib import Path

import numpy as np
from flint import arb
from numba import njit
from resonant_cycle import OUT, T, linear_trace
from scipy.integrate import solve_ivp


@njit
def field(t, q, gain):
    h = q[0] + 1j * q[1]
    b = 0.25 * np.sin(t) - 4 * np.sin(3 * t) + 3j * np.sin(6 * t)
    p = h * h * (3 + h) + b - 2
    dp = 3 * h * (2 + h)
    f = -gain * dp.conjugate() * p
    lin = -gain * abs(dp) ** 2
    anti = -6 * gain * (1 + h).conjugate() * p
    result = np.empty_like(q)
    result[0], result[1] = f.real, f.imag
    for i in (2, 4):
        v = q[i] + 1j * q[i + 1]
        d = lin * v + anti * v.conjugate()
        result[i], result[i + 1] = d.real, d.imag
    result[6] = 2 * lin
    return result


def orbit(gain, w, frame, tight):
    q0 = np.r_[frame @ w, frame[:, 0], frame[:, 1], 0.0]
    scale = np.linalg.norm(w)
    options = {
        "rtol": 2.3e-14 if tight else 3e-13,
        "atol": np.r_[
            np.repeat(scale * 2e-16, 2), np.full(4, 2e-14 if tight else 1e-13), 1e-21
        ],
        "max_step": 0.035 if tight else 0.09,
    }
    sol = solve_ivp(
        field,
        (0, 4 * T),
        q0,
        args=(gain,),
        method="DOP853",
        t_eval=np.arange(5) * T,
        **options,
    )
    if not sol.success:
        raise ArithmeticError(sol.message)
    phases = np.linalg.solve(frame, sol.y[:2]).T
    derivative = np.linalg.solve(frame, sol.y[2:6, -1].reshape(2, 2).T)
    return phases, derivative, float(sol.y[6, -1])


def spectrum(matrix, logdet):
    # Use the traceless part near the identity. Liouville fixes the determinant.
    disc = (matrix[0, 0] - matrix[1, 1]) ** 2 + 4 * matrix[0, 1] * matrix[1, 0]
    if disc < 0:
        mod = np.exp(logdet / 2)
        imag = np.sqrt(-disc) / 2
        real = np.sqrt(mod * mod - imag * imag)
        return "attracting", [[real, imag], [real, -imag]], mod
    hi = (np.trace(matrix) + np.sqrt(disc)) / 2
    lo = np.exp(logdet) / hi
    if not lo < 1 < hi:
        raise ArithmeticError("Unresolved saddle")
    return "saddle", [[lo, 0], [hi, 0]], hi


def experiment():
    row = json.loads((OUT / "resonances.json").read_text())["resonances"][-1]
    center = float(row["gain_center"])
    A, B = [
        complex(*(float(arb(v)) for v in row[k]))
        for k in ("cubic_isotropic", "cubic_resonant")
    ]
    sl, ht = (float(arb(v)) for v in row["frame_second_column"])
    frame = np.array([[1.0, sl], [0.0, ht]])
    g = np.sqrt(abs(B) ** 2 - A.real**2)
    records = []
    for offset in (2e-8, 8e-8, 3.2e-7):
        gain = center + offset
        delta = 4 * np.arccos(linear_trace(gain) / 2) - 2 * np.pi
        pair = []
        for sign in (-1, 1):
            theta = (np.angle(B) - np.arctan2(sign * g, -A.real)) / 4
            radius = np.sqrt(-delta / (A.imag + sign * g))
            w = radius * np.array([np.cos(theta), np.sin(theta)])
            for _ in range(10):
                phases, M, ld = orbit(gain, w, frame, False)
                residual = phases[-1] - w
                if np.linalg.norm(residual) < 3e-9 * radius:
                    break
                w -= np.linalg.solve(M - np.eye(2), residual)
            check, D, logdet = orbit(gain, w, frame, True)
            label, multipliers, mod = spectrum(D, logdet)
            oldlabel, oldvalues, _ = spectrum(M, ld)
            disagreement = np.max(np.abs(np.array(oldvalues) - multipliers))
            relclosure = np.linalg.norm(check[-1] - w) / radius
            radii = np.linalg.norm(check[:-1], axis=1) / radius
            sep = min(
                np.linalg.norm(a - b)
                for i, a in enumerate(check[:-1])
                for b in check[i + 1 : -1]
            )
            if not (
                oldlabel == label
                and disagreement < 1e-5
                and relclosure < 1e-5
                and np.min(radii) > 0.5
                and sep > radius
            ):
                raise ArithmeticError(
                    f"Unresolved cycle: {sign}, {offset}, {label}, {disagreement}, {relclosure}, {sep / radius}"
                )
            pair.append(
                {
                    "type": label,
                    "sign": sign,
                    "predicted_radius": radius,
                    "predicted_phase": theta,
                    "phases": check[:-1].tolist(),
                    "physical_offsets": (check[:-1] @ frame.T).tolist(),
                    "normalized_radii": radii.tolist(),
                    "relative_closure": relclosure,
                    "log_determinant": logdet,
                    "multipliers": multipliers,
                    "largest_multiplier_modulus": mod,
                    "solver_multiplier_disagreement": disagreement,
                }
            )
            print(
                offset,
                sign,
                label,
                "radius",
                radii,
                "res",
                relclosure,
                "logdet",
                logdet,
                "multipliers",
                multipliers,
                flush=True,
            )
        assert {x["type"] for x in pair} == {"saddle", "attracting"}
        records.append(
            {"gain": gain, "gain_offset": offset, "detuning": delta, "orbits": pair}
        )
    return {
        "case": 48000001,
        "resonance": row["gain_center"],
        "frame_second_column": [sl, ht],
        "records": records,
    }


def draw(data, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    paper, ink, muted = "#fcfbf7", "#16343e", "#516670"
    teal, red = "#087e78", "#b75a28"
    fig, axes = plt.subplots(
        1, 3, figsize=(9.0, 3.9), sharex=True, sharey=True, facecolor=paper
    )
    fig.subplots_adjust(left=0.09, right=0.98, top=0.69, bottom=0.22, wspace=0.14)
    fig.text(
        0.052,
        0.92,
        "Two four-cycles emerge beside an attracting point",
        size=20,
        weight="bold",
        color=ink,
    )
    fig.text(
        0.054,
        0.822,
        "Each dot is one return to the forcing phase. The same local coordinates and scale are used throughout.",
        size=10.5,
        color=muted,
    )
    for ax, record in zip(axes, data["records"]):
        ax.set_facecolor(paper)
        ax.axhline(0, lw=0.65, color="#d7dedc", zorder=0)
        ax.axvline(0, lw=0.65, color="#d7dedc", zorder=0)
        for orbit in record["orbits"]:
            points = np.array(orbit["phases"]) * 1e7
            closed = np.vstack([points, points[0]])
            sink = orbit["type"] == "attracting"
            color = teal if sink else red
            ax.plot(*closed.T, color=color, lw=1.4, ls="-" if sink else "--")
            ax.scatter(
                *points.T,
                marker="o" if sink else "D",
                s=28 if sink else 22,
                c=color if sink else paper,
                edgecolors=color,
                linewidth=1.25,
                zorder=3,
            )
        ax.scatter([0], [0], s=19, c=ink, marker="s", zorder=4)
        offset = record["gain_offset"] * 1e8
        ax.set_title(
            rf"$\Delta\kappa={offset:g}\times10^{{-8}}$", size=12, color=ink, pad=6
        )
        ax.set(
            aspect="equal",
            xlim=(-2.85, 2.85),
            ylim=(-2.85, 2.85),
            xticks=(-2, 0, 2),
            yticks=(-2, 0, 2),
        )
        ax.tick_params(labelsize=9, colors=muted, length=2)
        ax.spines[["top", "right"]].set_visible(False)
        for name in ("left", "bottom"):
            ax.spines[name].set_color("#bcc8c9")
    axes[0].set_ylabel(r"Im $w$ / $10^{-7}$", size=10, color=muted)
    axes[1].set_xlabel(r"Re $w$ / $10^{-7}$", size=10, color=muted, labelpad=1)
    handles = [
        Line2D([], [], c=ink, marker="s", ls="", ms=5, label="Stationary sink"),
        Line2D([], [], c=teal, marker="o", lw=1.5, ms=5, label="Attracting four-cycle"),
        Line2D(
            [],
            [],
            c=red,
            marker="D",
            mfc=paper,
            ls="--",
            lw=1.5,
            ms=5,
            label="Saddle four-cycle",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.008),
        ncol=3,
        frameon=False,
        fontsize=10,
        labelcolor=ink,
    )
    fig.savefig(path, dpi=175, facecolor=paper)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=OUT / "four-cycle.json")
    p.add_argument("--figure", type=Path, default=OUT / "four-cycle.png")
    args = p.parse_args()
    record = experiment()
    args.out.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    draw(record, args.figure)

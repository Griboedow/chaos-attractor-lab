"""Plot the measured saddle cycle and the certified leading drift at resonance."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PAPER, INK, MUTED = "#fbfaf6", "#203b43", "#52646d"
TEAL, RUST = "#087c77", "#ba512b"


def draw(data, path):
    rows = data["records"]
    coeff = complex(*data["resonant_coefficient"])
    vlim = (
        max(max(np.hypot(*v) for v in row["phase_zero_points"]) for row in rows) * 1.16
    )
    fig = plt.figure(figsize=(7.8, 5.3), facecolor=PAPER, layout=None)
    fig.text(
        0.055,
        0.912,
        "A saddle cycle passes through a sink",
        color=INK,
        fontsize=21,
        weight="bold",
    )
    fig.text(
        0.055,
        0.838,
        "Forced cubic  ·  middle window  ·  resonance at κ* ≈ 0.141705738258",
        color=MUTED,
        fontsize=11.5,
    )
    for i, side in enumerate((-1, 0, 1)):
        ax = fig.add_axes(
            [0.045 + 0.32 * i, 0.347, 0.274, 0.36], aspect="equal", facecolor=PAPER
        )
        ax.set(xlim=(-1.03, 1.03), ylim=(-1.03, 1.03))
        ax.set_axis_off()
        ax.set_title(
            ["Below κ*", "At κ*", "Above κ*"][i],
            size=13,
            weight="bold",
            color=INK,
            pad=4,
        )
        if side:
            for row in sorted(
                (r for r in rows if r["gain_offset"] * side > 0),
                key=lambda r: -abs(r["gain_offset"]),
            ):
                p = np.array(row["phase_zero_points"]) / vlim
                color = TEAL if side < 0 else RUST
                closed = np.vstack([p, p[0]])
                ax.plot(*closed.T, color=color, linewidth=0.95, alpha=0.50)
                ax.scatter(*p.T, s=23, color=color, zorder=3)
            ax.scatter(
                0, 0, marker="o", s=24, c=INK, edgecolors=PAPER, linewidth=0.6, zorder=4
            )
        else:
            phi = np.angle(coeff) / 3
            for j in range(6):
                v = np.exp(1j * (phi + j * np.pi / 3))
                ra, rb = (0.20, 0.90) if j % 2 == 0 else (0.90, 0.20)
                color = RUST if j % 2 == 0 else TEAL
                ax.annotate(
                    "",
                    xy=(rb * v.real, rb * v.imag),
                    xytext=(ra * v.real, ra * v.imag),
                    arrowprops={
                        "arrowstyle": "-|>",
                        "color": color,
                        "lw": 2,
                        "mutation_scale": 13,
                    },
                )
            ax.scatter(0, 0, s=24, facecolors=PAPER, edgecolors=INK, lw=1.4)
    fig.text(0.18, 0.291, "shrinking triangles", ha="center", color=TEAL, size=10.5)
    fig.text(0.50, 0.291, "three outward directions", ha="center", color=INK, size=10.5)
    fig.text(
        0.82, 0.291, "60° rotation in the limit", ha="center", color=RUST, size=10.5
    )
    fig.text(
        0.055,
        0.174,
        "Each triangle joins three consecutive returns to the same forcing phase.",
        color=INK,
        size=11.2,
    )
    fig.text(
        0.055,
        0.105,
        "The outer panels share the same linear coordinates and scale. The center arrows show",
        color=MUTED,
        size=10.2,
    )
    fig.text(
        0.055,
        0.057,
        "the leading three-period displacement at resonance.",
        color=MUTED,
        size=10.2,
    )
    fig.savefig(path, dpi=165, facecolor=PAPER)
    plt.close(fig)


if __name__ == "__main__":
    draw(json.loads((HERE / "three-cycle.json").read_text()), HERE / "three-cycle.png")

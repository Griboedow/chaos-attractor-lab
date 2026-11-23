"""Render three certified orbits of the same 18-by-18 matrix."""

import json

import matplotlib

matplotlib.use("Agg")
import itertools

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Rectangle
from sinkhorn_choreography import OUT


def matrix_phases(seed, sizes):
    labels = np.repeat(np.arange(3), sizes)
    n = len(labels)
    kernel = 0.99 * np.eye(n) + 0.01 / n
    x = np.array([0, *seed[:2]])[labels]
    y = np.array([0, *seed[2:]])[labels]
    table = np.exp(x[:, None] + y[None, :]) * kernel
    table /= table.sum()

    def step(matrix):
        matrix = matrix / matrix.sum(axis=1)[:, None] ** 1.75
        return matrix / matrix.sum(axis=0)[None, :] ** 1.75

    # Scaling the starting table by c scales its third return by c^((-3/4)^6).
    table *= step(step(step(table))).sum() ** (1 / (1 - (-0.75) ** 6))
    start = table.copy()
    for _ in range(3):
        yield table / table.sum(), labels
        table = step(table)
    np.testing.assert_allclose(table, start, atol=0, rtol=1e-11)


def main():
    packet = json.loads((OUT / "cycles.json").read_text())
    records = next(row for row in packet["dimensions"] if row["n"] == 18)
    by_size = {tuple(row["sizes"]): row for row in records["certificates"]}
    paper, ink, muted = "#faf9f5", "#162c3d", "#4f6271"
    accents = ["#21776c", "#c5522d", "#7060a5"]
    cmap = LinearSegmentedColormap.from_list(
        "mass", [paper, "#e5e8e4", "#99adb3", "#3b6272", "#112b41"]
    )
    norm = LogNorm(vmin=1e-5, vmax=1e-1, clip=True)
    fig = plt.figure(figsize=(8.5, 10.1), facecolor=paper)
    fig.text(
        0.07,
        0.962,
        "One matrix, three different orbits",
        size=25,
        weight="bold",
        color=ink,
    )
    fig.text(
        0.07,
        0.922,
        "Each 18 × 18 panel labels the share of mass in its dominant block.",
        size=13,
        color=muted,
    )
    lefts = [0.085, 0.387, 0.689]
    width = 0.22
    height = width * 8.5 / 10.1
    for i, left in enumerate(lefts):
        fig.text(
            left + width / 2,
            0.865,
            f"Step {i + 1}",
            size=15,
            weight="bold",
            ha="center",
            color=ink,
        )
    for sizes, bottom in zip([(6, 6, 6), (3, 6, 9), (1, 1, 16)], [0.602, 0.349, 0.096]):
        row = by_size[sizes]
        fig.text(
            0.07,
            bottom + height + 0.016,
            f"Groups of {sizes[0]}, {sizes[1]} and {sizes[2]} rows",
            size=13,
            color=ink,
            weight="bold",
        )
        for phase, ((matrix, groups), left) in enumerate(
            zip(matrix_phases(row["seed"], sizes), lefts)
        ):
            ax = fig.add_axes([left, bottom, width, height])
            ax.imshow(
                matrix, cmap=cmap, norm=norm, interpolation="nearest", origin="upper"
            )
            ax.set_axis_off()
            cuts = np.cumsum([0, *sizes])
            for j, (a, b) in enumerate(itertools.pairwise(cuts)):
                ax.add_patch(
                    Rectangle(
                        (a - 0.5, -1.17),
                        b - a,
                        0.49,
                        facecolor=accents[j],
                        edgecolor="none",
                        clip_on=False,
                    )
                )
                ax.add_patch(
                    Rectangle(
                        (-1.17, a - 0.5),
                        0.49,
                        b - a,
                        facecolor=accents[j],
                        edgecolor="none",
                        clip_on=False,
                    )
                )
            for cut in cuts[1:-1] - 0.5:
                ax.axvline(cut, color="#737f83", lw=0.6)
                ax.axhline(cut, color="#737f83", lw=0.6)
            ax.add_patch(
                Rectangle((-0.5, -0.5), 18, 18, fill=False, edgecolor="#879399", lw=0.7)
            )
            i, j = [(1, 0), (2, 1), (0, 2)][phase]
            fraction = matrix[np.ix_(groups == i, groups == j)].sum()
            fig.text(
                left + width / 2,
                bottom - 0.018,
                f"{fraction:.1%}",
                size=10.8,
                ha="center",
                color=muted,
            )
    ax = fig.add_axes([0.375, 0.031, 0.525, 0.011])
    ax.imshow(np.logspace(-5, -1, 256)[None, :], aspect="auto", cmap=cmap, norm=norm)
    ax.set_yticks([])
    ax.set_xticks([0, 64, 128, 192, 255], ["≤0.001%", "0.01%", "0.1%", "1%", "≥10%"])
    ax.tick_params(axis="x", labelsize=9, length=0, pad=4, colors=muted)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(0.07, 0.028, "Share in one cell", size=12, weight="bold", color=ink)
    fig.savefig(OUT / "cycles.png", dpi=140, facecolor=paper)
    plt.close(fig)


if __name__ == "__main__":
    main()

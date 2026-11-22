"""Render the observed motor itineraries on the cube of row orientations."""

import argparse
import json
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


def positions():
    even = np.array([[-1, -1, -1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]])
    vertices = np.stack([even, -even], axis=1).reshape((8, 3))
    horizontal = np.array([0.84, -0.54, 0.0])
    vertical = np.array([0.32, 0.50, 0.80])
    return vertices, vertices @ np.array([horizontal, vertical]).T


def draw_cube(ax, word, title, detail):
    vertices, xy = positions()
    for i, j in combinations(range(8), 2):
        if np.abs(vertices[i] - vertices[j]).sum() == 2:
            ax.plot(
                *xy[[i, j]].T, color="#acb8be", lw=1.6, solid_capstyle="round", zorder=0
            )
    ax.scatter(*xy.T, s=36, c="#a5b1b9", edgecolors="#ffffff", linewidths=0.6, zorder=1)
    route = list(word) + [word[0]]
    ax.plot(
        *xy[route].T,
        c="#087c72",
        lw=3.6,
        solid_capstyle="round",
        solid_joinstyle="round",
        zorder=2,
    )
    for index, pos in enumerate(word):
        x, y = xy[pos]
        ax.scatter([x], [y], s=230, c="#095f5c", zorder=3)
        ax.text(
            x,
            y,
            str(index + 1),
            ha="center",
            va="center",
            color="white",
            fontsize=9,
            zorder=4,
            weight="bold",
        )
    ax.set(aspect="equal", xlim=(-1.71, 1.71), ylim=(-1.85, 1.85))
    ax.axis("off")
    ax.text(
        0.5,
        1.02,
        title,
        transform=ax.transAxes,
        ha="center",
        fontsize=14,
        weight="bold",
        color="#172e39",
    )
    ax.text(
        0.5,
        -0.077,
        detail,
        transform=ax.transAxes,
        ha="center",
        fontsize=11,
        color="#40545f",
    )


def main(folder):
    runs = json.loads((folder / "historical.json").read_text())["runs"][1:]
    cycles = json.loads((folder / "cycles.json").read_text())
    fig = plt.figure(figsize=(10.2, 6.6), facecolor="white")
    fig.text(
        0.065,
        0.942,
        "Four motor rows, eight orientations",
        fontsize=24,
        color="#152f39",
        weight="bold",
    )
    fig.text(
        0.065,
        0.885,
        "Each vertex labels a copy of the thin braid. Numbers follow one periodic orbit.",
        fontsize=12,
        color="#40545f",
    )
    grid = fig.add_gridspec(
        1, 3, left=0.07, right=0.95, bottom=0.435, top=0.812, wspace=0.22
    )
    for j, (cycle, name) in enumerate(
        zip(cycles, ("One edge", "A triangle", "One square face"))
    ):
        ax = fig.add_subplot(grid[0, j])
        draw_cube(
            ax,
            cycle["checks"][1]["repeating_word"],
            name,
            f"$b={cycle['b']:g}$    $T={cycle['period']:.2f}$",
        )
    fig.add_artist(
        Line2D(
            [0.065, 0.95],
            [0.378, 0.378],
            color="#dbe1e3",
            lw=1,
            transform=fig.transFigure,
        )
    )
    fig.text(
        0.065,
        0.325,
        "Irregular orbits also favor cube edges",
        fontsize=16,
        weight="bold",
        color="#152f39",
    )
    fig.text(
        0.065,
        0.278,
        "Share of changes between successive classified visits",
        fontsize=11,
        color="#40545f",
    )
    ax = fig.add_axes([0.20, 0.055, 0.70, 0.18])
    colors = ("#087c72", "#e79c53", "#9daab4")
    for i, run in enumerate(runs):
        d = next(s for s in run["runs"] if "matrix" in s)
        left = 0
        for field, color in zip(("edge", "face", "opposite"), colors):
            width = d[field] / d["changes"]
            ax.barh(
                2 - i,
                width,
                height=0.53,
                left=left,
                color=color,
                edgecolor="white",
                linewidth=1,
            )
            if field == "edge":
                ax.text(
                    left + width * 0.50,
                    2 - i,
                    f"{width:.0%}",
                    ha="center",
                    va="center",
                    color="white",
                    weight="bold",
                    fontsize=10,
                )
            left += width
    ax.set(
        yticks=[2, 1, 0],
        yticklabels=[f"$b={r['b']:g}$" for r in runs],
        xlim=(0, 1),
        ylim=(-0.53, 2.5),
        xticks=[],
    )
    ax.tick_params(axis="y", length=0, pad=9, labelsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for x, color, label in zip(
        (0.200, 0.398, 0.674), colors, ("edge", "face diagonal", "body diagonal")
    ):
        fig.text(x, 0.026, "■", color=color, fontsize=11)
        fig.text(x + 0.019, 0.026, label, color="#40545f", fontsize=10)
    fig.savefig(folder / "motor-routes.png", dpi=170, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    main(p.parse_args().folder)

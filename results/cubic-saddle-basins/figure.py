"""Render phase-plane exits and a numerical slice through the three basins."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

COLORS = ("#087d80", "#bb4f2a", "#6952a0")


def sweep(w, turns, steps, keep=False):
    h = 2 * np.pi / steps
    t = np.linspace(0, 2 * np.pi, 2 * steps + 1)
    b = 0.25 * np.sin(t) - 4 * np.sin(3 * t) + 3j * np.sin(6 * t)
    path = []

    def field(v, drive):
        return -0.34 * (3 * v * (2 + v)).conj() * ((3 + v) * v * v - 2 + drive)

    for turn in range(turns):
        for j in range(steps):
            if keep and turn == turns - 1:
                path.append(1 + w.copy())
            first = field(w, b[2 * j])
            second = field(w + h / 2 * first, b[2 * j + 1])
            third = field(w + h / 2 * second, b[2 * j + 1])
            last = field(w + h * third, b[2 * j + 2])
            w = w + h / 6 * (first + 2 * second + 2 * third + last)
    if keep:
        path.append(1 + w.copy())
        return np.array(path)
    return 1 + w


def draws(proof, out):
    centers = np.array([complex(*p) for p in proof["sink_centers"]])
    unit = complex(*proof["unit_midpoint"])
    examples = proof["connections"]
    scale = np.array([float(r["scale"]) for r in examples])
    paths = sweep(scale * unit, 4, 2048, keep=True)
    target = np.array([r["target"] for r in examples])
    assert np.all(abs(paths[-1] - centers[target]) < proof["trapping_radius"])
    decades = np.linspace(4.0, 99.0, 1901)
    start = np.array([1, -1])[:, None] * unit * 10**-decades

    def classify(steps):
        endpoint = sweep(start.ravel(), 8, steps).reshape(start.shape)
        d = abs(endpoint[:, :, None] - centers[None, None, :])
        labels = d.argmin(axis=2)
        labels[d.min(axis=2) >= proof["trapping_radius"]] = 3
        return labels

    labels = classify(1024)
    coarse = classify(512)
    print(
        f"Grid: {labels.size} basin assignments; {np.count_nonzero(labels != coarse)} differ at half resolution."
    )
    fig = plt.figure(figsize=(9.6, 6.1), facecolor="white")
    fig.text(
        0.06,
        0.94,
        "Three basins meet at one saddle",
        fontsize=23,
        weight="bold",
        color="#172f3b",
    )
    fig.text(
        0.062,
        0.882,
        "Three tiny displacements from z = 1 reach different attracting cycles.",
        fontsize=12,
        color="#405563",
    )
    ax = fig.add_axes([0.09, 0.45, 0.82, 0.37])
    for idx, label in enumerate(target):
        curve = paths[:, idx]
        ax.plot(curve.real, curve.imag, lw=1.7, alpha=0.9, color=COLORS[label])
    for k, point in enumerate(centers):
        ax.scatter(
            point.real,
            point.imag,
            c=COLORS[k],
            s=43,
            edgecolor="white",
            linewidth=0.7,
            zorder=4,
        )
        shift = [(0.10, 0.09), (-0.10, 0.12), (0.05, -0.18)][k]
        ax.text(
            point.real + shift[0],
            point.imag + shift[1],
            "ABC"[k],
            c=COLORS[k],
            fontsize=12,
            weight="bold",
        )
    ax.scatter([1], [0], s=43, marker="D", color="#202f35", zorder=6)
    ax.annotate(
        "saddle",
        xy=(1, 0),
        xytext=(1.08, -0.34),
        color="#273d4a",
        size=10,
        arrowprops={"arrowstyle": "-", "color": "#63757e", "lw": 0.8},
    )
    ax.set(xlabel="Re z", ylabel="Im z", aspect="equal")
    ax.spines[["top", "right"]].set_visible(False)
    for name in ("left", "bottom"):
        ax.spines[name].set_color("#a4b1b7")
    ax.tick_params(labelsize=9, colors="#465b65")
    fig.text(
        0.062,
        0.335,
        "Moving closer to the saddle",
        fontsize=15,
        weight="bold",
        color="#172f3b",
    )
    fig.text(
        0.062,
        0.285,
        "Color gives the final cycle on either side of the unstable direction. Diamonds mark certified exits.",
        fontsize=10.5,
        color="#405563",
    )
    bands = fig.add_axes([0.12, 0.115, 0.83, 0.124])
    bands.imshow(
        labels,
        extent=(decades[0] - 0.025, decades[-1] + 0.025, -0.5, 1.5),
        interpolation="nearest",
        cmap=ListedColormap([*COLORS, "#b1bac0"]),
        vmin=0,
        vmax=3,
        aspect="auto",
        origin="upper",
    )
    for r in examples:
        s = float(r["scale"])
        bands.scatter(
            -np.log10(abs(s)),
            1 if s > 0 else 0,
            marker="D",
            s=38,
            c="white",
            edgecolor="#152e38",
            linewidth=0.9,
            zorder=3,
        )
    bands.axhline(0.5, color="white", lw=2)
    bands.set_yticks([0, 1], ["−", "+"])
    bands.set_xticks(
        [10, 20, 40, 60, 80, 95],
        [
            r"$10^{-10}$",
            r"$10^{-20}$",
            r"$10^{-40}$",
            r"$10^{-60}$",
            r"$10^{-80}$",
            r"$10^{-95}$",
        ],
    )
    bands.set_xlabel("Distance from z = 1", fontsize=10.5)
    bands.tick_params(length=0, pad=5, labelsize=10)
    for sp in bands.spines.values():
        sp.set_visible(False)
    fig.savefig(out, dpi=165, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("proof", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    draws(json.loads(args.proof.read_text()), args.out)

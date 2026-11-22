"""Render the frozen field-plane density and pulse timing from one native orbit."""

import argparse
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.signal import find_peaks


def main(a):
    z = np.load(a.cloud)["points"].astype("f8")
    peaks = find_peaks(z[:, 2], height=-1.50, distance=10, prominence=0.10)[0]
    ages = (
        np.arange(len(z))
        - peaks[
            np.maximum(np.searchsorted(peaks, np.arange(len(z)), side="right") - 1, 0)
        ]
    )
    k = np.arange(peaks[0], len(z))
    x, y, t = z[k].T
    lim = np.stack((x.min(), x.max(), y.min(), y.max()))
    xb = np.linspace(lim[0] - 0.3, lim[1] + 0.3, 1040)
    yb = np.linspace(lim[2] - 0.3, lim[3] + 0.3, 1220)
    h = np.histogram2d(x, y, bins=[xb, yb])[0].T
    av = np.histogram2d(x, y, bins=[xb, yb], weights=np.minimum(ages[k], 26))[
        0
    ].T / np.maximum(h, 1)
    bg = "#101923"
    base = np.array([16, 25, 35]) / 255
    ink = "#e8eef1"
    muted = "#a1b3c1"
    pal = LinearSegmentedColormap.from_list(
        "pulse age", ["#ffac6e", "#ffda99", "#cee6c9", "#79d7da", "#85a9f2"]
    )
    alpha = (1 - np.exp(-0.63 * h**0.57))[:, :, None]
    pixels = base * (1 - alpha) + pal(av / 26)[:, :, :3] * alpha
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "text.color": ink,
            "axes.labelcolor": muted,
            "ytick.color": muted,
            "xtick.color": muted,
            "axes.edgecolor": "#465466",
            "axes.facecolor": bg,
            "figure.facecolor": bg,
            "font.size": 10,
        }
    )
    fig = plt.figure(figsize=(12.4, 7.3), dpi=150, layout=None)
    fig.text(0.055, 0.934, "The scalloped cavity fan", fontsize=25, weight="bold")
    fig.text(
        0.056,
        0.891,
        "One orbit of a driven optical field with dynamic detuning",
        color=muted,
        fontsize=12,
    )
    ax = fig.add_axes([0.051, 0.099, 0.506, 0.737])
    ax.imshow(
        pixels,
        origin="lower",
        extent=[xb[0], xb[-1], yb[0], yb[-1]],
        interpolation="nearest",
        aspect="equal",
    )
    ax.set_xlabel("Real field amplitude")
    ax.set_ylabel("Imaginary field amplitude")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax = fig.add_axes([0.602, 0.571, 0.348, 0.213])
    start = 1200
    end = start + 107
    n = np.arange(start, end)
    hit = peaks[(peaks >= start) & (peaks < end)]
    ax.plot(n - start, z[n, 2], color="#91dce0", lw=1.15)
    ax.scatter(hit - start, z[hit, 2], s=22, color="#ffac6e", zorder=4)
    ax.set(xlim=(0, end - start - 1), xlabel="Iteration", ylabel="Detuning")
    ax.set_title("Detuning pulses", loc="left", fontsize=14, pad=8)
    ax = fig.add_axes([0.602, 0.221, 0.348, 0.196])
    lags = np.diff(peaks)
    bins = np.arange(lags.min(), lags.max() + 2) - 0.5
    freq, edges = np.histogram(lags, bins)
    xs = (edges[1:] + edges[:-1]) / 2
    ax.bar(
        xs,
        100 * freq / freq.sum(),
        width=0.83,
        color=["#ffac6e" if j == 20 else "#87b4c5" for j in xs],
    )
    ax.set(
        xlim=(lags.min() - 0.7, lags.max() + 0.7),
        xlabel="Iterations between successive peaks",
        ylabel="Share (%)",
    )
    ax.set_title("Pulse recurrence", loc="left", fontsize=14, pad=8)
    for ax in fig.axes[1:]:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
    cb = fig.add_axes([0.605, 0.100, 0.342, 0.020])
    cb.imshow(
        np.linspace(0, 1, 240)[None, :], aspect="auto", cmap=pal, extent=[0, 26, 0, 1]
    )
    cb.set_yticks([])
    cb.set_xticks([0, 5, 10, 15, 20, 26], ["0", "5", "10", "15", "20", "26+"])
    cb.tick_params(axis="x", length=2, labelsize=8)
    cb.set_xlabel("Field color tracks iterations since the preceding pulse", fontsize=9)
    fig.savefig(a.out, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--cloud", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    main(p.parse_args())

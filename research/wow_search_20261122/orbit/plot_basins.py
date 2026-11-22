"""Compare both long-lived trajectories at identical parameters and field scale."""

import argparse
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def main(a):
    z = [
        np.load(a.cache / f"{k}.npz")["points"]
        for k in ("case-86780", "case-86780-upper")
    ]
    xb = np.linspace(0.13, 1.35, 650)
    yb = np.linspace(-0.6, 0.6, 650)
    base = "#101923"
    ink = "#e7eff2"
    muted = "#a7bdc8"
    colors = ["#96dfdf", "#ffbe88"]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "text.color": ink,
            "axes.labelcolor": muted,
            "ytick.color": muted,
            "xtick.color": muted,
            "axes.edgecolor": "#455565",
            "axes.facecolor": base,
            "figure.facecolor": base,
            "font.size": 10,
        }
    )
    fig = plt.figure(figsize=(12, 6.8), dpi=145)
    fig.text(
        0.058,
        0.925,
        "Two detuning ranges at one pump setting",
        fontsize=24,
        weight="bold",
    )
    fig.text(
        0.06,
        0.870,
        "Each panel follows one initial condition under the same two-time-scale cavity map.",
        color=muted,
        fontsize=12,
    )
    for j, (pts, title) in enumerate(zip(z, ["Lower detuning", "Upper detuning"])):
        ax = fig.add_axes([0.075 + j * 0.298, 0.19, 0.275, 0.55])
        h, _, _ = np.histogram2d(pts[:, 0], pts[:, 1], bins=(xb, yb))
        c = LinearSegmentedColormap.from_list("branch", [base, colors[j]])
        ax.imshow(
            1 - np.exp(-0.48 * h.T**0.53),
            extent=[xb[0], xb[-1], yb[0], yb[-1]],
            origin="lower",
            interpolation="nearest",
            aspect="equal",
            cmap=c,
            vmin=0,
            vmax=1,
        )
        ax.set_title(title, fontsize=15, loc="left", color=colors[j])
        ax.set_xlabel("Real field amplitude")
        if j == 0:
            ax.set_ylabel("Imaginary field amplitude")
        else:
            ax.tick_params(labelleft=False)
        ax.spines[["top", "right"]].set_visible(False)
    ax = fig.add_axes([0.715, 0.32, 0.245, 0.34])
    for pts, col in zip(z, colors):
        ax.plot(np.arange(1400), pts[:1400, 2], lw=0.8, color=col)
    ax.set(xlim=(0, 1399), ylim=(4.85, 11.1), xlabel="Iteration", ylabel="Detuning")
    ax.axhline(7.5, ls="--", lw=0.65, color=muted, alpha=0.65)
    ax.set_title("Time traces", loc="left", fontsize=15)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.15)
    fig.text(
        0.06,
        0.115,
        "786,432 consecutive field samples per trajectory. Field-plane axes use identical absolute coordinates.",
        color=muted,
        fontsize=10,
    )
    fig.savefig(a.out, dpi=145)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    main(p.parse_args())

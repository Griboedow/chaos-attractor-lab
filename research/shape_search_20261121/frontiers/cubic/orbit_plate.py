"""Draw the one-turn sinks and the six-frame itineraries of the certified saddle orbits."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from gradient import field
from matplotlib.colors import LinearSegmentedColormap
from model import forcing
from scipy.integrate import solve_bvp

ROOT = Path(__file__).resolve().parents[4]
PACKET = ROOT / "results/cubic-forced-orbits"
SOURCE = ROOT / "results/shape_search_20261121/frontiers/cubic/certificates.json"
PAPER, INK, MUTED = "#faf9f4", "#172e3c", "#506575"
BLUE, GOLD = "#256393", "#cf702d"


def draw_paths(nodes, row, gain):
    """Use periodic collocation instead of propagating an unstable initial point."""
    mesh = np.linspace(0, 2 * np.pi, nodes.shape[1] + 1)
    smooth = np.linspace(0, 2 * np.pi, 1201)

    def rhs(t, y):
        deriv = field(y[0] + 1j * y[1], forcing(t, row), gain)
        return np.vstack((deriv.real, deriv.imag))

    paths = []
    worst = 0.0
    for orbit in nodes:
        values = np.r_[orbit, orbit[0]]
        solution = solve_bvp(
            rhs,
            lambda a, b: a - b,
            mesh,
            [values.real, values.imag],
            tol=2e-8,
            max_nodes=16000,
        )
        if not solution.success:
            raise ArithmeticError(solution.message)
        xy = solution.sol(mesh)
        worst = max(worst, np.max(abs(xy[0] + 1j * xy[1] - values)))
        xy = solution.sol(smooth)
        paths.append(xy[0] + 1j * xy[1])
    if worst > 2e-8:
        raise ArithmeticError(
            f"Drawn path missed its certified node proposals by {worst:g}"
        )
    print(f"Greatest path-to-proposal distance at the 64 phases: {worst:.3g}")
    return np.array(paths)


def render(paths, word_data, out):
    fig = plt.figure(figsize=(6.8, 8.75), facecolor=PAPER)
    fig.text(
        0.068,
        0.947,
        "21 orbits in one forced cubic",
        fontweight="bold",
        fontsize=23,
        color=INK,
    )
    fig.text(
        0.07, 0.908, "One forcing period  ·  gain κ = 0.34", fontsize=12, color=MUTED
    )
    fig.text(
        0.07,
        0.849,
        "Three attracting orbits",
        fontsize=15.5,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.07,
        0.817,
        "Paths in the complex plane. The dot marks the same starting phase.",
        fontsize=10,
        color=MUTED,
    )
    for i, path in enumerate(paths[:3]):
        ax = fig.add_axes([0.065 + i * 0.30, 0.668, 0.27, 0.139], facecolor=PAPER)
        ax.set_aspect("equal")
        ax.axhline(0, color="#deded5", linewidth=0.7)
        ax.axvline(0, color="#deded5", linewidth=0.7)
        ax.plot(path.real, path.imag, color="#455968", linewidth=1.55)
        ax.scatter(
            path[0].real,
            path[0].imag,
            color=INK,
            s=24,
            edgecolors=PAPER,
            linewidths=0.65,
            zorder=5,
        )
        ax.set(xlim=(-2.35, 2.35), ylim=(-1.25, 1.05))
        ax.set_axis_off()
    fig.text(
        0.07,
        0.635,
        "Eighteen saddle orbits",
        fontsize=15.5,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.07,
        0.603,
        "One row per orbit. Dots mark six certified signs of its real coordinate.",
        fontsize=10,
        color=MUTED,
    )
    fig.text(0.07, 0.571, "Left  (−)", fontsize=10.5, color=BLUE, fontweight="bold")
    fig.text(0.215, 0.571, "Right  (+)", fontsize=10.5, color=GOLD, fontweight="bold")
    fig.text(0.47, 0.571, "The final dot repeats the first.", fontsize=10, color=MUTED)
    colors = LinearSegmentedColormap.from_list("real_position", [BLUE, "#e8e6dc", GOLD])
    ax = fig.add_axes([0.207, 0.122, 0.732, 0.413], facecolor=PAPER)
    rows = sorted(
        word_data["saddle_words"], key=lambda v: (v["sign_changes"], v["word"])
    )
    frames = np.r_[np.array(word_data["phase_indices"]) / 64, 1.0]
    positions = []
    last_group = -1
    for row in rows:
        group = row["sign_changes"]
        line = (
            0 if not positions else positions[-1] + (1 if last_group == group else 1.8)
        )
        positions.append(line)
        last_group = group
        path = paths[row["orbit"]].real
        ax.imshow(
            path[None, :],
            extent=(0, 1, line + 0.34, line - 0.34),
            aspect="auto",
            interpolation="bilinear",
            cmap=colors,
            vmin=-1,
            vmax=1,
        )
        signs = row["word"] + row["word"][0]
        ax.scatter(
            frames,
            [line] * len(frames),
            s=13,
            c=[GOLD if s == "+" else BLUE for s in signs],
            edgecolors=PAPER,
            linewidths=0.65,
            clip_on=False,
            zorder=3,
        )
    for group in (0, 2, 4, 6):
        selected = [
            pos for pos, row in zip(positions, rows) if row["sign_changes"] == group
        ]
        midpoint = (selected[0] + selected[-1]) / 2
        ax.text(
            -0.034,
            midpoint - 0.10,
            f"{group} changes",
            transform=ax.get_yaxis_transform(),
            va="center",
            ha="right",
            fontsize=10,
            color=INK,
        )
    for j, midpoint in enumerate((frames[:-1] + frames[1:]) / 2):
        ax.text(
            midpoint,
            -1.05,
            "+ → −" if j % 2 == 0 else "− → +",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color=MUTED,
            clip_on=False,
        )
    ax.set_xlim(-0.006, 1.006)
    ax.set_ylim(positions[-1] + 0.81, -0.81)
    ax.set_yticks([])
    ax.set_xticks([0, 0.5, 1], labels=["0", "½", "1"])
    ax.set_xlabel("Forcing periods", fontsize=10, color=MUTED, labelpad=3)
    ax.tick_params(axis="x", colors=MUTED, length=2.5, labelsize=9)
    ax.spines[:].set_visible(False)
    fig.text(
        0.07,
        0.051,
        "Stay or switch as marked. These saddles realize all 18 closed patterns.",
        fontsize=10.2,
        color=INK,
    )
    fig.text(
        0.07,
        0.026,
        "With 0 or 4 sign changes, the expanding direction reverses after one turn.",
        fontsize=10.2,
        color=INK,
    )
    fig.savefig(out, dpi=180, facecolor=PAPER)
    plt.close(fig)


def main(folder):
    with np.load(folder / "shooting-nodes.npz") as packet:
        nodes = packet["nodes"]
    proof = json.loads((folder / "orbits.json").read_text())
    if [v["orbit"] for v in proof["orbits"] if v["type"] == "sink"] != [0, 1, 2]:
        raise ValueError("Unexpected attracting-orbit indices")
    row = next(
        r for r in json.loads(SOURCE.read_text())["rows"] if r["case"] == proof["case"]
    )
    gain = proof["gain"][0] / proof["gain"][1]
    paths = draw_paths(nodes, row, gain)
    render(
        paths,
        json.loads((folder / "itineraries.json").read_text()),
        folder / "orbits.png",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", type=Path, default=PACKET)
    main(parser.parse_args().folder)

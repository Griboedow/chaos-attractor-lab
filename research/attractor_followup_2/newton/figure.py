"""Draw the two perpendicular slices from the certified witness locations."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Wedge

OUT = Path(__file__).resolve().parents[3] / "results/attractor_followup_2/newton"


def main():
    data = json.loads((OUT / "certificate.json").read_text())
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "font.family": "DejaVu Sans"})
    fig, (geometry, chart) = plt.subplots(1, 2, figsize=(10.2, 5.45),
                                         gridspec_kw={"width_ratios": [0.91, 1.35]})
    fig.subplots_adjust(left=.04, right=.98, top=.86, bottom=.20, wspace=.34)
    fig.suptitle("Fifteen outcomes east and west, one north and south",
                 fontsize=16, fontweight="bold", x=.5, y=.978)
    geometry.set(xlim=(-1.26, 1.26), ylim=(-1.23, 1.23), aspect="equal")
    geometry.axis("off")
    geometry.add_patch(Wedge((0, 0), 1, 78, 102, color="#79a9cf"))
    geometry.add_patch(Wedge((0, 0), 1, 258, 282, color="#d7a06b"))
    for n in range(6):
        geometry.add_patch(Circle((0, 0), (5/6)**n, color="#73818b", fill=False,
                                  alpha=.52 if n < 4 else .3, linewidth=.8))
    geometry.plot([-1, 1], [0, 0], color="#1a343a", lw=1.2)
    geometry.scatter([0], [0], color="#111111", s=13)
    geometry.text(0, 1.07, "Root 2", ha="center", color="#285f87", weight="bold")
    geometry.text(0, -1.15, "Root 5", ha="center", color="#895629", weight="bold")
    geometry.text(-1.05, .08, "W", ha="right", weight="bold")
    geometry.text(1.05, .08, "E", ha="left", weight="bold")
    geometry.text(0, -.155, "pole", ha="center", fontsize=9)
    geometry.set_title("Two sectors, each 24° wide", pad=8, fontsize=12)
    labels = [f"Root {j}" for j in range(1, 7)] + [f"Phase {j}" for j in range(9)]
    for side, color, offset, mark in [("east", "#12636a", -.13, "o"),
                                      ("west", "#bf5e29", .13, "D")]:
        rows = data["horizontal"][side]
        x = [1000 / (1000 + row["j"]/60) for row in rows]
        y = [row["landing_label"] + offset for row in rows]
        chart.scatter(x, y, c=color, s=30, marker=mark, label=side.capitalize(), zorder=4)
    chart.set(yticks=np.arange(15), yticklabels=labels, ylim=(14.7, -.7),
              xlim=(5/6-.003, 1.003), xticks=[5/6, .875, .925, .975, 1.0],
              xticklabels=["5/6", ".875", ".925", ".975", "1"],
              xlabel="Distance from pole / band outer radius")
    chart.grid(axis="y", color="#e5e6e4", linewidth=.6)
    chart.axhline(5.5, color="#aab1ac", linewidth=.8)
    chart.spines["left"].set_visible(False)
    chart.tick_params(axis="y", length=0, labelsize=8.5)
    chart.set_title("One certified interval per outcome and ray", fontsize=12, pad=8)
    chart.legend(loc="lower right", bbox_to_anchor=(1., -0.255), ncol=2,
                 frameon=False, handletextpad=.05, columnspacing=1.25)
    fig.savefig(OUT / "newton-directional-basins.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()

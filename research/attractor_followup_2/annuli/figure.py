"""Plot the measured radii and certified attractor membership of two starts."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/attractor_followup_2/annuli"
case = json.loads((ROOT / "results/wow_search_20261122/symmetry/annuli.json").read_text())["13745"]
paths = json.loads((OUT / "capture.json").read_text())
_, lam, alpha, beta, gamma, omega, delta, e, h = map(float.fromhex, case["binary64_hex"])
ink, accent, band = "#263946", "#a83d1c", "#e2ebee"
fig, axes = plt.subplots(2, 1, figsize=(7.6, 4.6), sharex=True, sharey=True)
for index, ax in enumerate(axes):
    x, y = case["entrance_tests"][index]["start"]
    k = paths[index]["start_step"]
    z = complex(x, y)
    radii = []
    for _ in range(9):
        radii.append(abs(z))
        z = ((complex(lam, omega) + complex(alpha, delta)*abs(z)**2 + beta*(z**17).real)*z
             + gamma*z.conjugate()**16 + complex(e, h)*z.conjugate()**33)
    for lower, upper in case["annuli"]:
        ax.axhspan(float(lower), float(upper), color=band, zorder=0)
    ax.plot(range(9), radii, color="#7c8990", lw=1.2, zorder=1)
    ax.scatter(range(k), radii[:k], ec=ink, fc="white", s=36, zorder=2)
    ax.scatter(range(k, 9), radii[k:], color=accent, s=36, zorder=2)
    ax.set_title(f"({x}, {y})   ·   first entry at step {k}", loc="left", color=ink, fontsize=10)
    ax.set_ylabel("Radius")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#aeb8bd")
    ax.set_yticks([0.3, 0.5, 0.7])
    ax.tick_params(colors=ink, length=3)
axes[-1].set(xticks=range(9), xlabel="Ordinary map updates", xlim=(-0.28, 8.28), ylim=(0.28, 0.81))
fig.legend(handles=[Patch(color=band, label="Trapping annuli"),
                    Line2D([], [], marker="o", ls="", mec=ink, mfc="white", label="Outside attractor"),
                    Line2D([], [], marker="o", ls="", color=accent, label="Inside its interior")],
           ncol=3, frameon=False, loc="upper center", fontsize=9)
fig.subplots_adjust(left=.10, right=.98, bottom=.13, top=.83, hspace=.48)
fig.savefig(OUT / "first-entry.png", dpi=155, facecolor="white")

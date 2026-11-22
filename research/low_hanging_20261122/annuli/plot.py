"""Draw the inverse-branch box distances from the exact decimal witness."""
from decimal import Decimal
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter


packet = Path(__file__).resolve().parents[3] / "results/low_hanging_20261122/annuli"
witness = json.loads((packet / "witness.json").read_text())
points = [[Decimal(x), Decimal(y)] for x, y in witness["centers"]]
radii = [Decimal(r) / 10**9 for r in witness["half_width_in_billionths"]]
distance = [max(abs(x-points[0][0]), abs(y-points[0][1])) for x, y in points[1:]]
steps = list(range(1, len(points)))
fig, ax = plt.subplots(figsize=(7.2, 3.8), layout="constrained")
ax.set_yscale("log")
ax.axhspan(1.5e-9, float(radii[0]), color="#e2eeef")
ax.axhline(float(radii[0]), color="#42666c", lw=1.1, ls="--")
ax.errorbar(steps, [float(d) for d in distance], yerr=[float(r) for r in radii[1:]],
            color="#174c65", ecolor="#bd631f", capsize=4, fmt="o-", lw=1.5, ms=4)
ax.text(1.25, 1.5e-8, "Inside the starting square", color="#366067", fontsize=10)
ax.set(xlim=(0.7, 9.3), ylim=(1.5e-9, 6e-3), xticks=steps,
       xlabel="Backward steps of the 18-iterate return",
       ylabel="Distance from the initial center (max norm)")
ax.yaxis.set_major_locator(LogLocator(base=10, numticks=9))
ax.yaxis.set_minor_formatter(NullFormatter())
ax.grid(axis="y", color="#aebdc2", alpha=0.35)
ax.spines[["top","right"]].set_visible(False)
fig.savefig(packet / "inverse-branch.png", dpi=160)

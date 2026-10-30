"""Illustrate archived motor-cycle data and exact certified spectral brackets.

No trajectory is integrated and no proof decision is made by this renderer.
The orbit is a floating-point illustration; root intervals come from the
validated certificate. Plot coordinates use ordinary floating-point values.
"""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"docs/assets/motor-cycle-certificate.png")
    args = parser.parse_args()
    orbit_path = ROOT/"results/motor-rows/cluster-cycles/cycle-17.npz"
    certificate_path = ROOT/"results/motor-proof/continuation-recheck-20261030/target.json"
    orbit = np.load(orbit_path)
    certificate = json.loads(certificate_path.read_text())
    if certificate["status"] != "periodic_orbit_and_floquet_count_certified":
        raise ValueError("The referenced certificate does not certify the spectrum")
    intervals = certificate["spectrum"]["transverse"]["A"]["root_brackets"]
    if not all(row["passed"] for row in intervals):
        raise ValueError("A root bracket is not certified")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(12.6, 4.8), layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=(1.05, 1.2))
    ax = fig.add_subplot(gs[0], projection="3d")
    h, y, z = orbit["base"][:3]
    ax.plot(h, y, 1-z, color="#147d92", lw=2)
    ax.scatter([h[0]], [y[0]], [1-z[0]], color="#20324a", s=30)
    ax.set(xlabel="h", ylabel="Y", zlabel="C = 1 − Z",
           title="Four-plus-four periodic motion\nArchived approximate orbit; not a proof enclosure")
    ax.view_init(elev=23, azim=42)
    ax.set_box_aspect((1.25, 1, .85))
    ax = fig.add_subplot(gs[1])
    ax.axvspan(-1, 1, color="#e5f2ec", zorder=0)
    ax.axvline(-1, ls="--", color="#6c777d", lw=1)
    ax.axvline(1, ls="--", color="#6c777d", lw=1)
    for i, row in enumerate(intervals):
        lo, hi = map(float, row["bracket"])
        color = "#b94836" if i == 0 else "#147d92"
        ax.plot([lo, hi], [2-i]*2, color=color, lw=7, solid_capstyle="butt")
        ax.scatter([(lo+hi)/2], [2-i], color=color, s=28)
        ax.annotate(f"({row['bracket'][0]}, {row['bracket'][1]})", ((lo+hi)/2, 2-i),
                    xytext=(0, 14), textcoords="offset points", ha="center", color=color)
    ax.set(xlim=(-4.2, 1.3), ylim=(-.6, 2.65), xlabel="Real transverse Floquet multiplier",
           yticks=[0, 1, 2], yticklabels=["Root 3 × 6", "Root 2 × 6", "Root 1 × 6"],
           title="Certified open root intervals\nEach transverse root occurs six times")
    ax.text(.5, -.23, "Full 23D count: 6 unstable · 1 flow · 16 stable", transform=ax.transAxes,
            ha="center", fontsize=12, color="#20324a", fontweight="bold")
    ax.grid(axis="x", alpha=.2)
    fig.suptitle("A certified periodic saddle in the motor-row model", fontsize=15, fontweight="bold")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180, facecolor="white")
    plt.close(fig)
    manifest = dict(description="Visualization only; no new integration or proof computation",
                    inputs={str(p.relative_to(ROOT)).replace('\\', '/'):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (orbit_path, certificate_path)},
                    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    output_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest())
    args.output.with_suffix(".json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(str(args.output).encode("ascii", "backslashreplace").decode("ascii"))


if __name__ == "__main__":
    main()

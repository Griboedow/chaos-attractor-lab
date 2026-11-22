"""Locate the certified squares on the archived orbit density for Case 13745."""

from decimal import Decimal
from hashlib import sha256
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[3]
ORIGINAL = ROOT / "results/wow_search_20261122/symmetry"
STUDY = ROOT / "results/low_hanging_20261122/annuli"


def main():
    witness = json.loads((STUDY / "witness.json").read_bytes())
    case = witness["case"]
    record = json.loads((ORIGINAL / "numerics.json").read_bytes())["cases"][str(case)]
    source = ORIGINAL / f"cyclic-{case}.png"
    assert sha256(source.read_bytes()).hexdigest() == record["png_sha256"]
    # This is the same extent used by icons.hist for the archived PNG.
    extent = 1.035 * record["radius_range"][1]
    origin, *_, ninth = [[Decimal(value) for value in c] for c in witness["centers"]]
    nano = Decimal(10) ** 9
    cx, cy = [float((v - u) * nano) for u, v in zip(origin, ninth)]
    h0 = witness["half_width_in_billionths"][0]
    h9 = witness["half_width_in_billionths"][9]
    assert abs(cx) + h9 < h0 and abs(cy) + h9 < h0

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.edgecolor": "#5b6a72",
            "axes.labelcolor": "#283b47",
            "text.color": "#172c38",
            "xtick.color": "#455761",
            "ytick.color": "#455761",
        }
    )
    fig, (orbit, boxes) = plt.subplots(1, 2, figsize=(10.3, 5.2), layout="constrained")
    for ax in (orbit, boxes):
        ax.set_aspect("equal")
        ax.tick_params(length=3)
    orbit.set_title("Two million visits from one orbit", loc="left", fontsize=14, pad=11)
    orbit.imshow(plt.imread(source), origin="upper", extent=(-extent, extent, -extent, extent))
    orbit.set(
        xlabel="Real part of z",
        ylabel="Imaginary part of z",
        xticks=[-0.75, -0.5, 0, 0.5, 0.75],
        yticks=[-0.75, -0.5, 0, 0.5, 0.75],
    )
    xy = tuple(float(v) for v in origin)
    gold = "#ffd575"
    orbit.plot(*xy, "o", ms=10, mfc="none", mec=gold, mew=1.5)
    orbit.annotate(
        "Start of the proof",
        xy=xy,
        xytext=(-0.22, 0.08),
        color=gold,
        arrowprops={"arrowstyle": "-", "shrinkB": 5, "color": gold, "lw": 1.3},
        ha="left",
        va="center",
        fontsize=10,
    )

    boxes.set_title("Certified return inside the starting square", loc="left", fontsize=14, pad=11)
    limit = 1.18 * h0
    boxes.set(
        xlim=(-limit, limit),
        ylim=(-limit, limit),
        facecolor="#f6f8f9",
        xlabel="Real offset from center (10⁻⁹)",
        ylabel="Imaginary offset (10⁻⁹)",
        xticks=[-100, -50, 0, 50, 100],
        yticks=[-100, -50, 0, 50, 100],
    )
    boxes.set_axisbelow(True)
    boxes.grid(color="#dbe1e3", linewidth=0.65)
    boxes.add_patch(
        Rectangle((-h0, -h0), 2 * h0, 2 * h0, facecolor="none", edgecolor="#26485b", lw=2.1)
    )
    boxes.add_patch(
        Rectangle(
            (cx - h9, cy - h9), 2 * h9, 2 * h9, facecolor="#d2ece5", edgecolor="#167b6b", lw=1.8
        )
    )
    boxes.plot(0, 0, marker="+", ms=7, color="#26485b", mew=1.4)
    boxes.text(0.88 * h0, -0.87 * h0, "Starting square B₀", ha="right", color="#26485b")
    boxes.annotate(
        "Return box B₉",
        xy=(cx, cy + h9),
        xytext=(-52, 77),
        color="#10685b",
        arrowprops={"arrowstyle": "-", "color": "#167b6b", "lw": 1},
        ha="left",
        va="center",
    )
    fig.savefig(STUDY / "branches.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()

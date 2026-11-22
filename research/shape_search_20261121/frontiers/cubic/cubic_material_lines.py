"""Advect one material loop through three periods of the forced cubic oscillator."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from gradient import field
from matplotlib.colors import ListedColormap
from model import forcing
from numba import njit

REPO = Path(__file__).resolve().parents[4]
SOURCE = REPO / "results/shape_search_20261121/frontiers/cubic"
CASE = 48000001
GAINS = (0.03, 0.34)
CUTS = (-0.5, 0.85)
COLORS = ("#355b93", "#bfccce", "#ce6c39")
INK = "#182c37"


@njit(cache=True)
def advance(points, b, gain, first, stop, dt):
    result = np.empty_like(points)
    for j, z in enumerate(points):
        for i in range(first, stop):
            k1 = field(z, b[2 * i], gain)
            k2 = field(z + dt * k1 / 2, b[2 * i + 1], gain)
            k3 = field(z + dt * k2 / 2, b[2 * i + 1], gain)
            k4 = field(z + dt * k3, b[2 * i + 2], gain)
            z += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        result[j] = z
    return result


def evolve_block(points, b, gain, first, stop, dt, spacing):
    """Test the midpoint of each chord, then remove alternate redundant vertices."""
    midpoints = (points + np.roll(points, -1)) / 2
    after = advance(points, b, gain, first, stop, dt)
    centers = advance(midpoints, b, gain, first, stop, dt)
    following = np.roll(after, -1)
    split = (abs(centers - after) + abs(following - centers) > spacing) | (
        abs(centers - (after + following) / 2) > 0.002 * spacing
    )
    polygon = np.empty(len(after) + np.count_nonzero(split), dtype=complex)
    indices = np.arange(len(after)) + np.cumsum(split) - split
    polygon[indices] = after
    polygon[indices[split] + 1] = centers[split]
    if len(polygon) > 512:
        before, following = np.roll(polygon, 1), np.roll(polygon, -1)
        length = abs(before - polygon) + abs(following - polygon)
        remove = (length < spacing / 4) & (
            length - abs(following - before) < 0.001 * spacing
        )
        remove[1::2] = False
        remove[0] = False  # Preserve the starting marker for the unrolled view.
        polygon = polygon[~remove]
    return polygon


def sample(polygon):
    closed = np.r_[polygon, polygon[0]]
    arc = np.r_[0.0, np.cumsum(abs(np.diff(closed)))]
    if not np.isfinite(arc[-1]):
        raise ValueError("Nonfinite polygon")
    crossings = {}
    for cut in CUTS:
        for offset in (-0.05, 0.0, 0.05):
            right = closed.real > cut + offset
            crossings[f"{cut + offset:.2f}"] = int(
                np.count_nonzero(right[1:] != right[:-1])
            )
    info = {"length": float(arc[-1]), "crossings": crossings, "vertices": len(polygon)}
    return info, arc, closed


def ellipse(vertices):
    angle = np.arange(vertices) * 2 * np.pi / vertices
    return -0.75 - 0.10j + 1.50 * np.cos(angle) + 1.10j * np.sin(angle)


def fixed_tracers(row, steps=1024):
    drive = forcing(np.linspace(0, 2 * np.pi, 2 * steps + 1), row)
    polygon = ellipse(32768)
    records = [sample(polygon)[0]]
    for _ in range(3):
        polygon = advance(polygon, drive, GAINS[1], 0, steps, 2 * np.pi / steps)
        records.append(sample(polygon)[0])
    return {"gain": GAINS[1], "steps_per_turn": steps, "measurements": records}


def evolve(row, gain, steps, initial_vertices, blocks, spacing, draw, periods=3):
    dt = 2 * np.pi / steps
    drive = forcing(np.linspace(0, 2 * np.pi, 2 * steps + 1), row)
    polygon = ellipse(initial_vertices)
    records, snapshots = [], []
    for period in range(periods + 1):
        if period:
            for first in range(0, steps, steps // blocks):
                polygon = evolve_block(
                    polygon, drive, gain, first, first + steps // blocks, dt, spacing
                )
        info, arc, closed = sample(polygon)
        records.append(info)
        if draw and period <= 2:
            # Preserve equally spaced samples along the folded curve, including repeats.
            distance = np.linspace(0, arc[-1], 32769)
            snapshots.append(np.interp(distance, arc, closed))
        print(
            f"gain={gain:g} steps={steps} periods={period} length={arc[-1]:.5f}",
            flush=True,
        )
    return records, snapshots


def draw_figure(snapshots, data, out):
    background = "#fbfaf6"
    fine = data[-1]["runs"]
    fig = plt.figure(figsize=(9, 6.0), facecolor=background)
    fig.text(
        0.055,
        0.935,
        "How many strands fit in one silhouette?",
        size=22,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.057,
        0.88,
        "One loop carried by the cubic flow. Dashed lines mark two cuts in the phase plane.",
        size=10.5,
        color=INK,
    )
    panels = (
        (snapshots[0][0], "Initial loop"),
        (snapshots[0][2], "Two periods  ·  κ = 0.03"),
        (snapshots[1][2], "Two periods  ·  κ = 0.34"),
    )
    for i, (z, title) in enumerate(panels):
        ax = fig.add_axes([0.055 + i * 0.305, 0.53, 0.275, 0.26], facecolor=background)
        ax.set_xlim(-2.40, 2.0)
        ax.set_ylim(-1.40, 1.35)
        ax.set_aspect("equal")
        ax.set_axis_off()
        for cut, label in zip(CUTS, ("A", "B")):
            ax.axvline(cut, color="#829195", lw=0.95, ls=(0, (2, 3)))
            ax.text(
                cut,
                1.14,
                label,
                color=INK,
                va="center",
                ha="center",
                size=8.5,
                bbox={"facecolor": background, "edgecolor": "none", "pad": 1},
            )
        ax.plot(z.real, z.imag, color=INK, lw=1.45)
        ax.set_title(title, size=12, color=INK, pad=11, weight="semibold")

    fig.text(
        0.057,
        0.475,
        "Walk around the folded loop at κ = 0.34",
        size=15,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.057,
        0.431,
        "Color tracks position from left to right. Each stripe spans its share of the total arc length.",
        size=9.8,
        color=INK,
    )
    for j, (name, color) in enumerate(
        zip(("Left of A", "Between cuts", "Right of B"), COLORS)
    ):
        fig.text(0.057 + j * 0.164, 0.38, "■", size=11, color=color)
        fig.text(0.077 + j * 0.164, 0.381, name, size=8.8, color=INK)
    for period, top in ((1, 0.304), (2, 0.222)):
        z = snapshots[1][period]
        ax = fig.add_axes([0.18, top, 0.59, 0.047])
        ax.imshow(
            np.digitize(z.real, CUTS)[None, :],
            aspect="auto",
            interpolation="nearest",
            cmap=ListedColormap(COLORS),
            vmin=0,
            vmax=2,
        )
        ax.set_axis_off()
        fig.text(
            0.057,
            top + 0.011,
            f"{period} period" + ("s" if period == 2 else ""),
            size=11,
            color=INK,
        )
        fig.text(
            0.792,
            top + 0.011,
            f"Length  {fine[1][period]['length']:,.1f}",
            size=10,
            color=INK,
        )
    strong = fine[1][3]["crossings"]
    weak = fine[0][3]["crossings"]
    fig.text(
        0.057,
        0.128,
        f"After three periods the loop crosses A {strong['-0.50']:,} times "
        f"and B {strong['0.85']:,} times.",
        size=13,
        color=INK,
        weight="bold",
    )
    fig.text(
        0.057,
        0.076,
        f"At κ = 0.03, the corresponding counts are {weak['-0.50']} and {weak['0.85']}.",
        size=10,
        color=INK,
    )
    fig.savefig(out, dpi=170, facecolor=background)
    plt.close(fig)


def check_counts(data):
    for low, high in zip(data[0]["runs"], data[1]["runs"]):
        for a, b in zip(low, high):
            if a["crossings"] != b["crossings"]:
                raise ValueError("Crossing counts change under refinement")
    # The archived crossings are chronological; successive twists multiply on the left.
    generators = {-1: np.array([[1, -1], [0, 1]]), 2: np.array([[1, 0], [-1, 1]])}
    certificates = json.loads((SOURCE / "finite-gains.json").read_text())["records"]
    braid_actions = []
    for record, run in zip(certificates, data[1]["runs"], strict=True):
        action = np.eye(2, dtype=int)
        for crossing in record["dyadic_crossing_word"]:
            action = generators[crossing] @ action
        slope = np.array([1, 0])
        gate_counts = []
        for period in run:
            expected = (2 * np.abs(slope)).tolist()
            for cut, count in zip(CUTS, expected):
                if any(
                    period["crossings"][f"{cut + shift:.2f}"] != count
                    for shift in (-0.05, 0.0, 0.05)
                ):
                    raise ValueError(
                        "A measured cut count differs from the braid prediction"
                    )
            gate_counts.append(expected)
            slope = action @ slope
        braid_actions.append({"matrix": action.tolist(), "crossings": gate_counts})
    return braid_actions


def main(out):
    row = next(
        q
        for q in json.loads((SOURCE / "certificates.json").read_text())["rows"]
        if q["case"] == CASE
    )
    resolutions = ((1024, 256, 128, 0.035), (2048, 512, 256, 0.0175))
    data, snapshots = [], []
    for steps, initial, blocks, spacing in resolutions:
        runs = []
        for gain in GAINS:
            records, views = evolve(
                row, gain, steps, initial, blocks, spacing, steps == 2048
            )
            runs.append(records)
            if views:
                snapshots.append(views)
        data.append(
            {
                "steps_per_turn": steps,
                "initial_vertices": initial,
                "remeshes_per_turn": blocks,
                "target_spacing": spacing,
                "runs": runs,
            }
        )
    braid_actions = check_counts(data)
    out.mkdir(parents=True, exist_ok=True)
    packet = {
        "case": CASE,
        "gains": list(GAINS),
        "center": [-0.75, -0.1],
        "semiaxes": [1.5, 1.1],
        "cuts": list(CUTS),
        "offset_checks": [-0.05, 0.0, 0.05],
        "chronological_braid_actions": braid_actions,
        "fixed_tracers": fixed_tracers(row),
        "resolutions": data,
    }
    (out / "measurements.json").write_text(
        json.dumps(packet, indent=2) + "\n", encoding="utf-8"
    )
    draw_figure(snapshots, data, out / "folds.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=REPO / "results/cubic-material-lines"
    )
    main(parser.parse_args().out)

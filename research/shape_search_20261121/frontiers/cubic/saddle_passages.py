"""Validate six phases of covering relations for the driven cubic using its boundary arcs."""

import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from functools import cache, lru_cache
from pathlib import Path

from flint import acb, acb_poly, arb, ctx
from periodic_boxes import coefficient_norm, conjugate, drive_series, local_polynomials

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "results/cubic-saddle-passages"
KNOTS = (0, 11, 21, 32, 43, 53, 64)  # units of 2*pi/64
# Columns give two local coordinate axes. Entries are exact multiples of 1/1000.
# Sites are ordered (-1, +1). Frames three phases later are their exact reflections.
FRAMES = (
    ((983, -519, -184, 855), (855, 184, 519, 983)),
    ((935, -208, -355, 978), (989, 656, 145, 755)),
    ((978, -355, -208, 935), (755, 145, 656, 989)),
)
SUBSTEPS = 12
DEGREE = 9
BITS = 80


@cache
def axes(phase, site):
    if phase >= 3:
        u, v = axes(phase - 3, -site)
        return u.conjugate(), -v.conjugate()
    a, b, c, d = FRAMES[phase][(site + 1) // 2]
    return acb(arb(a) / 1000, arb(c) / 1000), acb(arb(b) / 1000, arb(d) / 1000)


@cache
def forcing_cells(phase):
    h = 2 * arb.pi() / (64 * SUBSTEPS)
    return h, [
        drive_series(i * h, h, DEGREE)
        for i in range(KNOTS[phase] * SUBSTEPS, KNOTS[phase + 1] * SUBSTEPS)
    ]


def transport_cell(center, radius, b, step, gain):
    """Enclose the image of a Euclidean disk by a Taylor center and a log-norm error."""
    forcing, tail = b
    path, _ = local_polynomials(center, [], forcing, step, gain)
    derivative = 3 * path * path - 3
    residual = path * path * path - 3 * path + acb_poly(forcing)
    defect = coefficient_norm(
        path.derivative() + step * gain * conjugate(derivative) * residual
    )
    defect += step * gain * coefficient_norm(derivative) * tail

    # Test an a priori error radius by first exit over this time cell.
    cap = max(arb(1) / 1024, 2 * radius)
    span = coefficient_norm(path - acb_poly([center])) + cap
    size = abs(center) + span
    dp = abs(3 * center * center - 3)
    dp_min = max(arb(0), (dp - 6 * abs(center) * span - 3 * span**2).lower())
    p_max = (
        abs(center**3 - 3 * center + forcing[0])
        + dp * span
        + 3 * abs(center) * span**2
        + span**3
    )
    p_max += step * 121 / 4  # |b'| <= 1/4 + 12 + 18
    growth = (step * gain * (-(dp_min**2) + 6 * size * p_max)).exp().upper()
    peak = max(arb(1), growth)
    if not peak * (radius + defect) < cap:
        return None
    endpoint = path(1)
    exact_midpoint = endpoint.mid()
    radius = (growth * radius + peak * defect + abs(endpoint - exact_midpoint)).upper()
    if not radius < arb(3) / 20:
        return None
    return exact_midpoint, radius


def target_margin(center, radius, phase, target, edge, orientation, width):
    u, v = axes(phase, target)
    det = (u.conjugate() * v).imag
    if not det > 0:
        raise ArithmeticError("Invalid frame")
    z = center - target
    x = (z.conjugate() * v).imag / det + arb(0, radius * abs(v) / det)
    y = (u.conjugate() * z).imag / det + arb(0, radius * abs(u) / det)
    if edge < 2:
        return ((1 if edge else -1) * orientation * x - width).lower()
    # Avoid the two complete vertical rays beyond the target entry edges.
    return max((abs(x) - width).lower(), (width - abs(y)).lower())


def check_edge(job):
    phase, site, edge = job
    ctx.prec = BITS
    gain, width = arb(17) / 50, arb(7) / 50
    step, forcing = forcing_cells(phase)
    u, v = axes(phase, site)
    fixed, moving = (u, v) if edge < 2 else (v, u)
    base = site + (1 if edge % 2 else -1) * width * fixed
    switching_site = 1 if phase % 2 == 0 else -1
    targets = [site, -site] if site == switching_site else [site]
    absent = [] if site == switching_site else [-site]
    orientation = 1 if site == switching_site else -1
    pending = [(0, 0)]
    leaves, deepest = 0, 0
    margin = arb(10)
    while pending:
        level, index = pending.pop()
        param = (arb(2 * index + 1) / 2**level - 1) * width
        z = base + param * moving
        radius = (abs(moving) * width / 2**level + abs(z - z.mid())).upper()
        disk = (z.mid(), radius)
        for b in forcing:
            disk = transport_cell(*disk, b, step, gain)
            if disk is None:
                break
        if disk is not None:
            clearance = min(
                target_margin(*disk, phase + 1, t, edge, orientation, width)
                for t in targets
            )
            for t in absent:
                # A real-coordinate half-plane is convex. If the complete boundary
                # image lies in it, so does the embedded disk.
                z, rad = disk
                a, b = axes(phase + 1, t)
                separation = (
                    site * (z.real - t) - rad - width * (abs(a.real) + abs(b.real))
                )
                clearance = min(clearance, separation.lower())
            if clearance > 0:
                leaves += 1
                deepest = max(deepest, level)
                margin = min(margin, clearance)
                continue
        if level == 36:
            raise ArithmeticError(
                f"Unresolved edge {job} at dyadic interval {index}/2^{level}"
            )
        pending.extend(((level + 1, 2 * index + 1), (level + 1, 2 * index)))
    print(
        f"phase={phase} site={site:+d} edge={edge} leaves={leaves} depth={deepest}",
        flush=True,
    )
    return {
        "phase": phase,
        "site": site,
        "edge": edge,
        "targets": targets,
        "excluded": absent,
        "orientation": orientation,
        "dyadic_intervals": leaves,
        "maximum_depth": deepest,
        "minimum_clearance": str(margin),
    }


def sink_separation():
    import numpy as np

    source = ROOT / "results/cubic-forced-orbits"
    proof = json.loads((source / "orbits.json").read_text())
    path = source / "shooting-nodes.npz"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if (
        digest != proof["node_file_sha256"]
        or proof["box_half_width"] != "2^-40"
        or proof["case"] != 48000001
        or proof["gain"] != [17, 50]
    ):
        raise ValueError("Unexpected orbit certificate or node boxes")
    with np.load(path) as file:
        nodes = file["nodes"]
    distances = []
    for record in proof["orbits"]:
        if record["type"] == "sink":
            x = arb(float(nodes[record["orbit"], 0].real), arb(2) ** -40)
            margins = []
            for site in (-1, 1):
                u, v = axes(0, site)
                extent = (arb(7) / 50) * (abs(u.real) + abs(v.real))
                margins.append((abs(x - site) - extent).lower())
            distances.append(min(margins))
    if len(distances) != 3 or not all(gap > arb(3) / 10 for gap in distances):
        raise ArithmeticError("A sink was not separated from the two phase-zero boxes")
    return {
        "node_sha256": digest,
        "phase_zero_real_distance_bounds": [str(d) for d in distances],
    }


def route_counts(edges):
    import numpy as np

    steps = []
    for phase in range(3):
        matrix = np.zeros((2, 2), dtype=int)
        for site in (-1, 1):
            sides = [e for e in edges if (e["phase"], e["site"]) == (phase, site)]
            if sorted(e["edge"] for e in sides) != [0, 1, 2, 3]:
                raise ValueError("Each claimed passage requires all four source sides")
            targets = set.intersection(*(set(e["targets"]) for e in sides))
            absent = set.intersection(*(set(e["excluded"]) for e in sides))
            if targets & absent or targets | absent != {-1, 1}:
                raise ValueError("Unclassified destination")
            for target in targets:
                matrix[(site + 1) // 2, (target + 1) // 2] = 1
        steps.append(matrix)
    product = np.eye(2, dtype=int)
    for matrix in steps + [m[::-1, ::-1] for m in steps]:
        product = product @ matrix
    first, second = int(product.trace()), int((product @ product).trace())
    if product.tolist() != [[5, 8], [8, 13]]:
        raise ValueError("Unexpected passage counts")
    return {
        "one_period": product.tolist(),
        "closed_words": [first, second],
        "least_two_period_orbits_lower": (second - first) // 2,
    }


def figure(packet_path):
    import matplotlib.pyplot as plt
    import numpy as np
    from cubic_material_lines import evolve_block
    from model import forcing

    packet = json.loads(packet_path.read_text())
    if len(packet["edges"]) != 24:
        raise ValueError("Incomplete edge packet")
    source = ROOT / "results/shape_search_20261121/frontiers/cubic/certificates.json"
    row = next(
        r for r in json.loads(source.read_text())["rows"] if r["case"] == 48000001
    )
    steps = 4096
    drive = forcing(np.linspace(0, 2 * np.pi, 2 * steps + 1), row)
    paper, ink, soft = "#faf9f4", "#172e3c", "#506575"
    colors = {-1: "#24658f", 1: "#bd602a"}
    fig = plt.figure(figsize=(8.8, 5.6), facecolor=paper)
    fig.text(
        0.055, 0.925, "Passages between two saddles", size=22, weight="bold", color=ink
    )
    fig.text(
        0.057,
        0.865,
        "First interval of the forcing cycle. Each colored curve outlines an evolved box.",
        size=11,
        color=soft,
    )

    def box(phase, site):
        u, v = (complex(float(c.real), float(c.imag)) for c in axes(phase, site))
        return site + 0.14 * np.array([-u - v, u - v, u + v, -u + v])

    for index, site in enumerate((-1, 1)):
        ax = fig.add_axes([0.055 + 0.47 * index, 0.55, 0.42, 0.20], facecolor=paper)
        center = box(0, site)
        z = np.concatenate(
            [
                np.linspace(a, b, 64, endpoint=False)
                for a, b in zip(center, np.roll(center, -1))
            ]
        )
        for first in range(0, KNOTS[1] * steps // 64, 4):
            z = evolve_block(z, drive, 0.34, first, first + 4, 2 * np.pi / steps, 0.006)
        for target in (-1, 1):
            corners = box(1, target)
            ax.fill(
                corners.real,
                corners.imag,
                facecolor="#ecece3",
                edgecolor=ink,
                linewidth=1.2,
            )
            ax.text(target, 0.30, f"{target:+d}", ha="center", color=soft, size=10)
        ax.plot(z.real, z.imag, color=colors[site], lw=2)
        ax.set(xlim=(-1.87, 1.88), ylim=(-0.59, 0.53), aspect="equal")
        ax.set_axis_off()
        ax.set_title(
            f"Start in the {site:+d} box".replace("-", "−"),
            loc="left",
            fontsize=13,
            color=colors[site],
            fontweight="bold",
        )
    fig.text(
        0.057,
        0.445,
        "Allowed passages over one forcing period",
        size=15,
        weight="bold",
        color=ink,
    )
    ax = fig.add_axes([0.12, 0.225, 0.825, 0.16], facecolor=paper)
    times = np.asarray(KNOTS) / 64
    for site in (-1, 1):
        ax.plot(times, np.full(7, site), color="#bbc4c4", lw=1.5)
        ax.scatter(times, np.full(7, site), color=colors[site], s=42, zorder=4)
        ax.text(
            -0.03,
            site,
            f"{site:+d}".replace("-", "−"),
            size=11,
            color=colors[site],
            va="center",
            ha="right",
            weight="bold",
        )
    for j in range(6):
        site = 1 if j % 2 == 0 else -1
        ax.annotate(
            "",
            xy=(times[j + 1], -site),
            xytext=(times[j], site),
            arrowprops={
                "arrowstyle": "-|>",
                "mutation_scale": 13,
                "color": colors[site],
                "lw": 1.5,
                "shrinkA": 4,
                "shrinkB": 5,
            },
        )
    ax.set(
        xlim=(-0.02, 1.02),
        ylim=(-1.45, 1.45),
        yticks=[],
        xticks=[0, 0.5, 1],
        xticklabels=["0", "½ period", "1 period"],
    )
    ax.tick_params(length=0, labelcolor=soft)
    ax.spines[:].set_visible(False)
    fig.text(
        0.057,
        0.093,
        "Every allowed itinerary occurs for some initial point, however long the route.",
        fontsize=12.5,
        color=ink,
        weight="bold",
    )
    fig.savefig(packet_path.with_name("passages.png"), dpi=165, facecolor=paper)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--out", type=Path, default=OUT / "coverings.json")
    parser.add_argument(
        "--figure-only", action="store_true", help="draw from an existing certificate"
    )
    args = parser.parse_args()
    ctx.prec = BITS
    if args.figure_only:
        figure(args.out)
        return
    jobs = [(p, site, edge) for p in range(3) for site in (-1, 1) for edge in range(4)]
    if args.workers == 1:
        records = list(map(check_edge, jobs))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            records = list(executor.map(check_edge, jobs))
    packet = {
        "gain": "17/50",
        "phases_over_64": KNOTS,
        "frames_times_1000": FRAMES,
        "half_width": "7/50",
        "degree": DEGREE,
        "substeps_per_64th": SUBSTEPS,
        "arb_bits": BITS,
        "edges": records,
        "known_sinks": sink_separation(),
        "route_counts": route_counts(records),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(packet, indent=2) + "\n", encoding="utf-8")
    figure(args.out)


if __name__ == "__main__":
    main()

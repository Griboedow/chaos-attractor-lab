"""Draw the phase hexagons in native amplitude coordinates and their unfolded incidence."""

import argparse
from collections import deque
from pathlib import Path
import math
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Circle, Polygon
import numpy as np
from verify import compute

PAPER = "#f4efe5"
INK = "#182d38"
FOIL = "#f3ca79"
CORNERS = "#324a55"
HEX = "#397f7a"
PENT1 = "#367672"
PENT2 = "#664d70"
EDGE = "#c4e3c9"
CHORD = "#ffca77"
METRIC = np.diag([-1.0, 1.0, 1.0])


def normalize(p):
    return p / np.sqrt(-np.sum((p @ METRIC) * p, axis=-1, keepdims=True))


def disk(p):
    radius = np.linalg.norm(p[..., 1:], axis=-1, keepdims=True)

    return p[..., 1:] * np.tanh(np.arcsinh(radius) / 3) / np.maximum(radius, 1e-15)


def angle_points(t, r, angles):
    return np.column_stack(
        (np.full(len(angles), t), r * np.cos(angles), r * np.sin(angles))
    )


def geodesic(a, b, n=16):
    u = np.linspace(0, 1, n)[:, None]
    return normalize((1 - u) * a + u * b)


def boundary(corners):
    return np.vstack(
        [geodesic(a, b)[:-1] for a, b in zip(corners, np.roll(corners, -1, axis=0))]
    )


def phase_tile():
    angles = np.arange(6) * np.pi / 3 + np.deg2rad(21)
    v = angle_points(np.sqrt(3), np.sqrt(2), angles)
    m = angle_points(np.sqrt(2), 1.0, angles + np.pi / 6)
    p = normalize(m + np.roll(m, 1, axis=0))
    normals = angle_points(1.0, np.sqrt(2), angles + np.pi / 6)
    mirrors = [np.eye(3) - 2 * np.outer(n, n @ METRIC) for n in normals]
    return v, m, p, mirrors


def chambers(mirrors, limit=6.55):
    queue = deque([(np.eye(3), 1)])
    seen = set()
    while queue:
        g, parity = queue.popleft()
        center = g[:, 0]
        key = tuple(np.rint(center[1:] * 1e7).astype(int))
        if key in seen or center[0] > np.cosh(limit):
            continue
        seen.add(key)
        yield g, parity
        for h in mirrors:
            queue.append((g @ h, -parity))


def add_disk(ax, equal):
    v, m, p, mirrors = phase_tile()
    pops = [(boundary(m), HEX)]
    if equal:
        pops = [
            (boundary([p[0], m[0], m[1], m[2], p[3]]), PENT1),
            (boundary([p[3], m[3], m[4], m[5], p[0]]), PENT2),
        ]
    path_inner = [geodesic(m[i - 1], m[i]) for i in range(6)]
    diameter = geodesic(v[0], v[3], 28)
    fills = []
    colors = []
    edges = []
    cuts = []
    seams = []
    for g, parity in chambers(mirrors):
        for poly, col in pops:
            fills.append(disk(poly @ g.T))
            colors.append(col)
        for path in path_inner:
            edges.append(disk(path @ g.T))
        if equal:
            cuts.append(disk(diameter @ g.T))
        if parity > 0:
            for i in range(6):
                seams.append(disk(geodesic(v[i], v[(i + 1) % 6]) @ g.T))
    radius = 0.941
    clip = Circle((0, 0), radius, facecolor=CORNERS, edgecolor="none")
    ax.add_patch(clip)

    def layer(collection):
        collection.set_clip_path(clip)
        ax.add_collection(collection)

    layer(PolyCollection(fills, facecolors=colors, edgecolors="none", zorder=2))
    layer(
        LineCollection(seams, colors="#6b7e80", linewidths=0.37, alpha=0.64, zorder=3)
    )
    layer(LineCollection(edges, colors=EDGE, linewidths=0.88, zorder=4))
    if equal:
        layer(LineCollection(cuts, colors=CHORD, linewidths=1.04, zorder=5))
    ax.add_patch(Circle((0, 0), radius, fill=False, edgecolor=INK, lw=1.5, zorder=6))
    ax.set(xlim=(-0.958, 0.958), ylim=(-0.958, 0.958), aspect="equal")
    ax.axis("off")


def physical(ax, b):
    c = 3 / 8
    p = (c * c - 1 - b * b) / 2
    h = math.sqrt(b * b - p * p)
    # Actual corners are intersections of successive amplitude bounds.
    normals = [(1, 0), (-1, -1), (0, 1), (1, 0), (-1, -1), (0, 1)]
    rhs = [1, -c, b, -1, c, -b]
    corners = [
        np.linalg.solve([normals[i - 1], normals[i]], [rhs[i - 1], rhs[i]])
        for i in range(6)
    ]
    corners = [(x, -x - y) for x, y in corners]
    ax.add_patch(
        Polygon(corners, closed=True, facecolor="#ebe2cf", edgecolor="#567078", lw=1)
    )
    t = np.linspace(0, 2 * np.pi, 600)
    curve = np.array([np.cos(t), p * np.cos(t) - h * np.sin(t)])
    ax.plot(curve[0], -curve[0] - curve[1], color="#287a71", lw=1.8)
    if b == 1:
        ax.plot([-1, 1], [0, 0], color="#b67221", lw=1.7)
    for n in [(1, 0), (0, 1), (-1, -1)]:
        u = np.array([[1, p], [p, b * b]]) @ n
        u /= math.sqrt(n @ u)
        ax.scatter(
            [u[0], -u[0]], [-u.sum(), u.sum()], s=9, color=INK, zorder=4, linewidth=0
        )
    ax.set(xlim=(-1.18, 1.18), ylim=(-0.52, 0.52), aspect="equal")
    ax.axis("off")


def figure(path):
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "text.color": INK, "figure.facecolor": PAPER}
    )
    fig = plt.figure(figsize=(12.8, 8.3), dpi=170)
    fig.text(
        0.061,
        0.928,
        "One diameter doubles a geodesic mosaic",
        fontsize=25,
        weight="bold",
    )
    fig.text(
        0.063,
        0.885,
        "Straight rays on a periodic surface of genus fifteen",
        fontsize=13,
        color="#47616a",
    )
    for j, (b, title) in enumerate([(9 / 8, "Unequal weights"), (1, "Equal weights")]):
        center = 0.269 + 0.466 * j
        fig.text(
            center,
            0.828,
            title + ("   b = 1" if b == 1 else "   b = 9/8"),
            ha="center",
            fontsize=15,
            weight="bold",
        )
        ax = fig.add_axes([center - 0.213, 0.224, 0.426, 0.592])
        add_disk(ax, b == 1)
        physical(fig.add_axes([center - 0.208, 0.107, 0.119, 0.100]), b)
        res = compute(7, b == 1)
        counts = res["faces_by_edges"]
        typ = 5 if b == 1 else 6
        txt = (
            f"{res['rays']} rays  ·  {res['vertices']} crossings\n{counts[4]} quadrilaterals + {counts[typ]} "
            + ("pentagons" if b == 1 else "hexagons")
        )
        fig.text(
            center - 0.075, 0.156, txt, fontsize=11.5, va="center", linespacing=1.5
        )
    fig.text(
        0.5,
        0.048,
        "Pale arcs lift the ellipse. Gold lifts its diameter. Insets show the cosine amplitudes A and C at equal scale.",
        ha="center",
        fontsize=10,
        color="#526872",
    )
    fig.savefig(path, dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    figure(ap.parse_args().out)

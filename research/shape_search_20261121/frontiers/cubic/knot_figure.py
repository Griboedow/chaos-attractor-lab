"""Draw the measured period-four attractor in an untwisted periodic Floquet frame."""

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
from four_cycle import OUT, T, field
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LinearSegmentedColormap
from scipy.integrate import solve_ivp

PAPER, INK, MUTED = "#fbfaf5", "#153b3f", "#4d6263"
COLORS = LinearSegmentedColormap.from_list(
    "lap", ["#087d83", "#388eaa", "#ba6b48", "#c49447", "#087d83"], N=1024
)


def orbit_in_tube(record, n=1024):
    gain = record["gain"]
    orbit = next(o for o in record["orbits"] if o["type"] == "attracting")
    phase = np.linspace(0, T, n + 1)
    t = np.linspace(0, 4 * T, 4 * n + 1)

    def solve(q, times):
        sol = solve_ivp(
            field,
            (0, times[-1]),
            q,
            args=(gain,),
            t_eval=times,
            method="DOP853",
            rtol=2.4e-14,
            atol=[2e-22, 2e-22, 2e-15, 2e-15, 2e-15, 2e-15, 1e-21],
            max_step=0.04,
        )
        if not sol.success:
            raise ArithmeticError(sol.message)
        return sol.y

    q = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0])
    fundamental = solve(q, phase)[2:6].T.reshape(-1, 2, 2).transpose(0, 2, 1)
    monodromy = fundamental[-1]
    theta = np.arccos(np.trace(monodromy) / 2)
    e1 = np.array([1.0, 0.0])
    basis = np.column_stack(
        [e1, (np.cos(theta) * e1 - monodromy[:, 0]) / np.sin(theta)]
    )
    # X(T) C = C R(-theta). Lift that rotation to the physical clockwise turn
    # observed during a forcing period, and undo it to make a periodic frame.
    a = phase * (1 + theta / T)
    rot = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]]).transpose(2, 0, 1)
    frame = fundamental @ basis @ rot
    assert np.linalg.det(basis) > 0 and np.allclose(frame[0], frame[-1], atol=2e-6)
    # Polar rotation angle of an orientation-preserving real 2-by-2 matrix.
    polar = frame[:, 0, 0] + frame[:, 1, 1] + 1j * (frame[:, 1, 0] - frame[:, 0, 1])
    frame_turns = np.unwrap(np.angle(polar)) / T
    assert np.max(abs(np.diff(frame_turns))) < 0.05
    assert abs(frame_turns[-1] - frame_turns[0]) < 1e-6

    q[:2] = orbit["physical_offsets"][0]
    h = solve(q, t)[:2].T
    w = np.linalg.solve(frame[np.arange(len(t)) % n], h[..., None])[:, :, 0]
    w = (w[:, 0] + 1j * w[:, 1]) / orbit["predicted_radius"]
    physical_angle = np.unwrap(np.angle(h[:, 0] + 1j * h[:, 1]))
    turns = (physical_angle - physical_angle[0]) / T
    model = np.exp(1j * (np.angle(w[0]) - 1.25 * t))
    # Four equally spaced model strands have separation sqrt(2) in each fiber.
    deviation = np.max(abs(w - model))
    assert abs(turns[-1] + 5) < 1e-6 and deviation < 0.2
    print(
        f"physical turns {turns[-1]:.12f}; normalized model discrepancy {deviation:.6g}"
    )
    return t, w, turns


def embedded(t, w):
    radial = 2.30 + w.real
    return np.column_stack([radial * np.cos(t), radial * np.sin(t), w.imag])


def camera(azimuth, elevation):
    az, el = np.radians([azimuth, elevation])
    eye = np.array([np.cos(az) * np.cos(el), np.sin(az) * np.cos(el), np.sin(el)])
    right = np.array([-np.sin(az), np.cos(az), 0])
    return np.array([right, np.cross(eye, right), eye])


def catalog_inputs(record, folder):
    t, w, _ = orbit_in_tube(record, n=2048)
    curve = embedded(t, w)
    standard = camera(-45, 27)
    poses = {"overhead": [20, 74], "oblique": [20, 16]}
    views = {}
    for key, angles in poses.items():
        rotation = camera(*angles).T @ standard
        assert (
            np.allclose(rotation @ rotation.T, np.eye(3))
            and np.linalg.det(rotation) > 0
        )
        views[key] = (curve @ rotation).astype("<f4")
    folder.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(folder / "knot-views.npz", **views)
    meta = {
        "description": "Two rigid views of one four-period orbit. All 8193 uniform-time states include both endpoints.",
        "embedding": "((2.3 + Re w) cos t, (2.3 + Re w) sin t, Im w) in the normalized periodic Floquet frame",
        "gain": record["gain"],
        "start_azimuth_elevation_degrees": poses,
        "center": [0, 0, 0],
        "basis": np.eye(3).tolist(),
        "projection": standard[:2].tolist(),
    }
    (folder / "knot-views.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8", newline="\n"
    )


def draw_knot(ax, t, w):
    curve = embedded(t, w)[:-1:2]
    tangent = np.roll(curve, -1, axis=0) - np.roll(curve, 1, axis=0)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    normal = np.cross(tangent, [0, 0, 1])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    binormal = np.cross(tangent, normal)
    angle = 2 * np.pi * np.arange(20) / 20
    normals = (
        normal[:, None] * np.cos(angle)[None, :, None]
        + binormal[:, None] * np.sin(angle)[None, :, None]
    )
    mesh = curve[:, None] + 0.069 * normals

    right, up, eye = camera(20, 74)
    projected = mesh @ np.column_stack([right, up, eye])
    corners = np.stack(
        [
            projected,
            np.roll(projected, -1, 0),
            np.roll(projected, (-1, -1), (0, 1)),
            np.roll(projected, -1, 1),
        ],
        axis=2,
    ).reshape(-1, 4, 3)
    depth = corners[:, :, 2].mean(axis=1)
    order = np.argsort(depth)
    light = -0.40 * right + 0.68 * up + 0.61 * eye
    light /= np.linalg.norm(light)
    half = (eye + light) / np.linalg.norm(eye + light)
    faces = normals + np.roll(normals, -1, axis=1)
    faces /= np.linalg.norm(faces, axis=-1)[..., None]
    rgb = COLORS(t[:-1:2] / (4 * T))[:, None, :3]
    diffuse = 0.58 + 0.43 * np.maximum(faces @ light, 0)
    glint = 0.31 * np.maximum(faces @ half, 0) ** 16
    rgb = np.clip(rgb * diffuse[..., None] + glint[..., None], 0, 1).reshape(-1, 3)
    ax.add_collection(
        PolyCollection(
            corners[order, :, :2],
            facecolors=rgb[order],
            edgecolors="face",
            linewidths=0.19,
            antialiased=True,
        )
    )
    x, y = corners[..., 0], corners[..., 1]
    ax.set(
        xlim=(x.min() - 0.10, x.max() + 0.10),
        ylim=(y.min() - 0.10, y.max() + 0.10),
        aspect="equal",
    )
    ax.axis("off")


def draw(record, path):
    t, w, turns = orbit_in_tube(record)
    fig = plt.figure(figsize=(11.0, 5.55), facecolor=PAPER)
    fig.text(
        0.045,
        0.908,
        "Four forcing periods, five clockwise turns",
        size=25,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.047,
        0.845,
        "The attracting four-cycle of the forced cubic near its third 1:4 resonance",
        size=13,
        color=MUTED,
    )
    ax = fig.add_axes([0.017, 0.105, 0.602, 0.706], facecolor=PAPER)
    draw_knot(ax, t, w)
    fig.text(
        0.316,
        0.057,
        "The (4, −5) torus knot in a periodic Floquet frame",
        ha="center",
        size=10.5,
        color=MUTED,
    )

    chart = fig.add_axes([0.663, 0.256, 0.292, 0.442], facecolor=PAPER)
    chart.set_title(
        "Turns about $z=1$", loc="left", color=INK, size=15, weight="bold", pad=13
    )
    p = np.column_stack([t / T, turns])
    segments = np.stack([p[:-1], p[1:]], axis=1)
    chart.add_collection(
        LineCollection(
            segments, colors=COLORS(t[:-1] / (4 * T)), linewidths=2.9, capstyle="round"
        )
    )
    returns = p[:: (len(p) - 1) // 4]
    chart.scatter(
        *returns.T,
        c=COLORS(returns[:, 0] / 4),
        s=40,
        edgecolors=PAPER,
        linewidths=1,
        zorder=3,
        clip_on=False,
    )
    chart.set(
        xlim=(-0.04, 4.04),
        ylim=(-5.2, 0.2),
        xticks=range(5),
        yticks=range(-5, 1),
        xlabel="Elapsed forcing periods",
    )
    chart.xaxis.label.set(color=MUTED, size=10.5)
    chart.xaxis.labelpad = 7
    chart.tick_params(colors=MUTED, labelsize=10, length=0, pad=6)
    chart.grid(axis="both", color="#dbe0d9", lw=0.65, zorder=0)
    chart.set_axisbelow(True)
    for spine in chart.spines.values():
        spine.set_visible(False)
    fig.text(
        0.663,
        0.102,
        "Angle measured directly in the physical plane",
        size=10,
        color=MUTED,
    )
    fig.text(0.663, 0.057, rf"$\kappa = {record['gain']:.13f}$", size=10, color=MUTED)
    fig.savefig(path, dpi=180, facecolor=PAPER)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", type=Path, default=OUT / "four-cycle.json")
    output = parser.add_mutually_exclusive_group()
    output.add_argument("--out", type=Path, default=OUT / "torus-knot.png")
    output.add_argument("--catalog-inputs", type=Path, metavar="DIRECTORY")
    args = parser.parse_args()
    record = json.loads(args.packet.read_text())["records"][0]
    if args.catalog_inputs:
        catalog_inputs(record, args.catalog_inputs)
    else:
        draw(record, args.out)

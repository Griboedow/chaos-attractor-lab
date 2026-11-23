"""Map the cubic basins in Cartesian and log-polar coordinates; keep pixel grids outside Git."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
PACKET = ROOT / "results/cubic-saddle-basins"
COLORS = ["#087d80", "#bb4f2a", "#6952a0", "#f6f4ea"]
EXTENT = [-2.25, 2.25, -1.44, 1.44]
RADIUS_LIMITS = [2.0, 80.0]  # -log10(r)
PERIODS = 8
ENTRY_RADIUS = 0.004

KERNEL = r"""
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double2 rhs(double2 w, double2 forcing, double site) {
    double2 square = (double2)(w.x*w.x-w.y*w.y, 2*w.x*w.y);
    double2 deriv = 6*site*w + 3*square;
    double2 cubic = (double2)(square.x*w.x-square.y*w.y, square.x*w.y+square.y*w.x);
    double2 poly = forcing + (double2)(-2*site, 0.) + 3*site*square + cubic;
    return -.34*(double2)(deriv.x*poly.x+deriv.y*poly.y, deriv.x*poly.y-deriv.y*poly.x);
}
__kernel void classify(__global const double2 *offsets, __global const char *anchors,
                       __global const double2 *drive, __global const double2 *sinks,
                       __global uchar *destinations, __global uchar *turns,
                       int steps, int max_turns, double h, double radius) {
    int idx = get_global_id(0);
    double2 w = offsets[idx];
    double s = (double)anchors[idx];
    uchar label = 3;
    int turn;
    for (turn=1; turn<=max_turns; turn++) {
        for (int j=0; j<steps; j++) {
            double2 a = rhs(w, drive[2*j], s);
            double2 b = rhs(w+h*.5*a, drive[2*j+1], s);
            double2 c = rhs(w+h*.5*b, drive[2*j+1], s);
            double2 d = rhs(w+h*c, drive[2*j+2], s);
            w += (h/6)*(a+2*b+2*c+d);
            double next = (w.x+s >= 0) ? 1. : -1.;
            w.x += s-next;
            s = next;
        }
        for (int k=0; k<3; k++) {
            double2 diff = w-(sinks[k]-(double2)(s,0.));
            if (dot(diff,diff) < radius*radius) { label=k; break; }
        }
        if (label!=3 || !isfinite(w.x) || !isfinite(w.y)) break;
    }
    destinations[idx] = label;
    turns[idx] = label==3 ? 0 : turn;
}
"""


def drive(t):
    return 0.25 * np.sin(t) - 4 * np.sin(3 * t) + 3j * np.sin(6 * t)


def centers():
    return json.loads((PACKET / "certificate.json").read_text())["sink_centers"]


class Device:
    def __init__(self):
        import pyopencl as cl

        devices = [
            d
            for p in cl.get_platforms()
            for d in p.get_devices(device_type=cl.device_type.GPU)
            if "cl_khr_fp64" in d.extensions
        ]
        self.cl = cl
        self.ctx = cl.Context([devices[0]])
        self.q = cl.CommandQueue(self.ctx)
        self.program = cl.Program(self.ctx, KERNEL).build()
        self.kernel = cl.Kernel(self.program, "classify")
        self.sinks = self.buf(np.array(centers(), dtype=np.float64))
        self.name = devices[0].name

    def buf(self, a):
        flags = self.cl.mem_flags
        return self.cl.Buffer(
            self.ctx,
            flags.READ_ONLY | flags.COPY_HOST_PTR,
            hostbuf=np.ascontiguousarray(a),
        )

    def integrate(self, w, anchors, steps):
        cl = self.cl
        b = self.buf(drive(np.arange(2 * steps + 1) * np.pi / steps))
        labels = np.empty(w.size, np.uint8)
        turns = np.empty_like(labels)
        offsets, sites = w.ravel(), anchors.ravel()
        for i in range(0, w.size, 8192):
            chunk = slice(i, i + 8192)
            size = offsets[chunk].size
            lab = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size)
            cap = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size)
            self.kernel(
                self.q,
                (size,),
                None,
                self.buf(offsets[chunk]),
                self.buf(sites[chunk]),
                b,
                self.sinks,
                lab,
                cap,
                np.int32(steps),
                np.int32(PERIODS),
                np.float64(2 * np.pi / steps),
                np.float64(ENTRY_RADIUS),
            )
            cl.enqueue_copy(self.q, labels[chunk], lab)
            cl.enqueue_copy(self.q, turns[chunk], cap)
        return labels.reshape(w.shape), turns.reshape(w.shape)

    def checked(self, w, s):
        coarse, _ = self.integrate(w, s, 1024)
        fine, capture = self.integrate(w, s, 2048)
        differences = coarse != fine
        refinements = []
        steps = 4096
        while differences.any() and steps <= 16384:
            new, new_turn = self.integrate(w[differences], s[differences], steps)
            indices = np.flatnonzero(differences)
            disagree = new != fine.flat[indices]
            refinements.append([steps, len(indices), int(disagree.sum())])
            fine.flat[indices] = new
            capture.flat[indices] = new_turn
            differences.flat[indices] = disagree
            steps *= 2
        fine[differences] = 3
        return fine, {
            "samples": fine.size,
            "refinements": refinements,
            "counts": np.bincount(fine.ravel(), minlength=4).tolist(),
            "latest_capture": int(capture.max()),
        }


def grid():
    x = EXTENT[0] + (np.arange(800) + 0.5) * (EXTENT[1] - EXTENT[0]) / 800
    y = EXTENT[2] + (np.arange(512) + 0.5) * (EXTENT[3] - EXTENT[2]) / 512
    z = x[None, :] + 1j * y[:, None]
    s = np.where(z.real >= 0, 1, -1).astype(np.int8)
    theta = (np.arange(720) + 0.5) * 2 * np.pi / 720
    d = 2 + (np.arange(624) + 0.5) / 8
    polar = 10.0 ** (-d[:, None]) * np.exp(1j * theta[None, :])
    return z - s, s, polar, d


def directions():
    # The reverse-time stationary tangent is the small real eigenvector of M.
    # Compute it as the annihilator of the dominant left eigenvector to avoid
    # resolving the reciprocal sub-unit multiplier of this ill-conditioned matrix.
    n = 16384
    h = 2 * np.pi / n
    b = drive(np.arange(2 * n + 1) * h / 2)
    conj_coeff = 12 * 0.34 - 6 * 0.34 * b
    a = np.array(
        [[conj_coeff.real, conj_coeff.imag], [conj_coeff.imag, -conj_coeff.real]]
    ).transpose(2, 0, 1)
    m = np.eye(2)
    for j in range(n):
        k1 = a[2 * j] @ m
        k2 = a[2 * j + 1] @ (m + h * k1 / 2)
        k3 = a[2 * j + 1] @ (m + h * k2 / 2)
        k4 = a[2 * j + 2] @ (m + h * k3)
        m += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    eigenvalues, eigenvectors = np.linalg.eig(m.T)
    i = np.argmax(abs(eigenvalues))
    fixed = json.loads((PACKET / "certificate.json").read_text())
    lam = float(fixed["eigenvalue"].split()[0].lstrip("["))
    if abs(eigenvalues[i] / lam - 1) > 1e-8:
        raise ArithmeticError(
            "The independent linear return misses the certified expanding multiplier"
        )
    normal = complex(*eigenvectors[:, i].real)
    v = 1j * normal
    if v.imag < 0:
        v = -v
    u = complex(*fixed["unit_midpoint"])
    return u, v, lam


def cross(a, b):
    return (np.conjugate(a) * b).imag


def compute(cache):
    gpu = Device()
    start, sites, polar, decades = grid()
    ones = np.ones(polar.shape, np.int8)
    atlas, global_stats = gpu.checked(start, sites)
    local, polar_stats = gpu.checked(polar, ones)
    u, v, lam = directions()
    shorter = decades < 60
    shifted, shifted_stats = gpu.checked(polar[shorter] / lam, ones[shorter])
    projected_w = u * cross(polar, v) / cross(u, v)
    projected, projected_stats = gpu.checked(projected_w, ones)
    # Concentrate samples next to the stable tangent. Deviations of order r^2
    # deliberately probe where the tangent approximation loses its accuracy.
    q = -1 + (np.arange(160) + 0.5) * 2 / 160
    r = 10.0 ** (-(2 + (np.arange(80) + 0.5) * 4 / 80))[:, None]
    transverse = u * q[None, :] * r**2
    stress = np.array([sign * v * r + transverse for sign in (-1, 1)])
    full_near, _ = gpu.checked(stress, np.ones(stress.shape, np.int8))
    pred_near, _ = gpu.checked(transverse, np.ones(transverse.shape, np.int8))
    stats = {
        "gain": "17/50",
        "phase": 0,
        "classification_radius": ENTRY_RADIUS,
        "trapping_radius": 0.0055,
        "device": gpu.name,
        "cartesian_extent": EXTENT,
        "cartesian": global_stats,
        "polar": polar_stats,
        "theta_degrees": [0, 360],
        "negative_log10_radius": RADIUS_LIMITS,
        "stable_tangent_degrees": float(np.degrees(np.angle(v))),
        "polar_grid_shape": list(polar.shape),
        "cartesian_grid_shape": list(start.shape),
        "unstable_direction": [u.real, u.imag],
        "stable_direction": [v.real, v.imag],
        "expanding_multiplier": lam,
        "scaled_check": dict(
            domain_decades=[2, 60],
            **shifted_stats,
            different_labels=int((shifted != local[shorter]).sum()),
        ),
        "projection_check": dict(
            **projected_stats, different_labels=int((projected != local).sum())
        ),
        "close_to_tangent": {
            "samples": stress.size,
            "different_projection_labels": int(
                np.count_nonzero(full_near != pred_near[None, :, :])
            ),
            "unresolved": int(
                np.count_nonzero((full_near == 3) | (pred_near[None, :, :] == 3))
            ),
            "starts": "w=±r*v+q*r^2*u. v stable tangent, u unstable, q in (-1,1), 2 < -log10(r) < 6.",
        },
    }
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, atlas=atlas, local=local, stats=json.dumps(stats))
    print(json.dumps(stats), flush=True)


def cpu_check(w, s):
    from scipy.integrate import solve_ivp

    capture = np.full(w.shape, 3, dtype=np.uint8)
    z, anchors = w.copy(), s.copy()
    targets = np.array([complex(*q) for q in centers()])
    for _ in range(PERIODS):
        todo = np.flatnonzero(capture == 3)
        if not len(todo):
            break
        y, sites = z[todo], anchors[todo]
        for j in range(32):

            def field(t, h, sites=sites):
                return (
                    -0.34
                    * (6 * sites * h + 3 * h**2).conj()
                    * (drive(t) - 2 * sites + 3 * sites * h**2 + h**3)
                )

            answer = solve_ivp(
                field,
                (j * np.pi / 16, (j + 1) * np.pi / 16),
                y,
                method="DOP853",
                rtol=2e-12,
                atol=np.maximum(1e-140, np.minimum(abs(y), 1) * 2e-14),
            )
            if not answer.success:
                raise ArithmeticError(answer.message)
            y = answer.y[:, -1]
            new = np.where(y.real + sites >= 0, 1, -1).astype(np.int8)
            y += sites - new
            sites = new
        z[todo], anchors[todo] = y, sites
        distance = abs(y[:, None] - (targets[None, :] - sites[:, None]))
        keep = distance.min(axis=1) < ENTRY_RADIUS
        capture[todo[keep]] = distance.argmin(axis=1)[keep]
    return capture


def render(cache):
    from matplotlib import pyplot as plt
    from matplotlib.colors import ListedColormap

    with np.load(cache) as file:
        atlas, local = file["atlas"], file["local"]
        stats = json.loads(str(file["stats"]))
    w, s, polar, _ = grid()
    rng = np.random.default_rng(48000001)
    trials = []
    for start, sites, labels in (
        (w, s, atlas),
        (polar, np.ones(polar.shape, np.int8), local),
    ):
        transition = np.zeros(labels.shape, bool)
        transition[1:, :] |= labels[1:] != labels[:-1]
        transition[:, 1:] |= labels[:, 1:] != labels[:, :-1]
        ix = np.unique(
            np.r_[
                rng.choice(start.size, 128, replace=False),
                rng.choice(np.flatnonzero(transition), 128, replace=False),
            ]
        )
        obtained = cpu_check(start.flat[ix], sites.flat[ix])
        if np.any(obtained != labels.flat[ix]):
            raise ArithmeticError(
                "Independent DOP853 sample differs from the displayed destinations"
            )
        trials.append({"samples": len(ix), "mismatches": 0})
    stats["cpu_dop853_checks"] = trials
    stats["cpu_sample_order"] = ["cartesian", "polar"]
    stats["count_order"] = ["A", "B", "C", "unresolved"]
    stats["refinement_columns"] = [
        "steps per turn",
        "retested",
        "disagreeing with previous resolution",
    ]
    stats["pixel_order"] = (
        "Rows start at the lower imaginary endpoint (Cartesian) or largest radius (polar)."
    )
    for key, data in (("cartesian", atlas), ("polar", local)):
        stats[key]["pixel_sha256"] = hashlib.sha256(data.tobytes()).hexdigest()
    (PACKET / "basin-map.json").write_text(
        json.dumps(stats, indent=2) + "\n", encoding="utf-8"
    )

    paper, ink = "#fbfaf6", "#233943"
    fig = plt.figure(figsize=(6.6, 8.0), facecolor=paper)
    fig.text(
        0.065,
        0.943,
        "Three basins meet at a saddle",
        size=21.5,
        weight="bold",
        color=ink,
    )
    fig.text(
        0.067,
        0.907,
        "Initial states at phase zero  ·  forced cubic, κ = 0.34",
        size=11,
        color="#526675",
    )
    fig.text(
        0.067,
        0.875,
        "Color gives the final cycle A, B or C. Diamonds mark the two saddles.",
        size=9.8,
        color="#526675",
    )
    cmap = ListedColormap(COLORS)
    ax = fig.add_axes([0.12, 0.541, 0.79, 0.318])
    ax.imshow(
        atlas,
        extent=EXTENT,
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=3,
    )
    ax.set(xlabel="Re z", ylabel="Im z", xticks=[-2, -1, 0, 1, 2], yticks=[-1, 0, 1])
    for i, (x, y) in enumerate(centers()):
        ax.scatter(x, y, s=165, facecolor=COLORS[i], edgecolor=paper, linewidth=1.3)
        ax.text(
            x,
            y,
            "ABC"[i],
            ha="center",
            va="center",
            size=8.4,
            color="white",
            weight="bold",
        )
    for x in (-1, 1):
        ax.scatter(
            x, 0, s=42, marker="D", facecolor=ink, edgecolor=paper, linewidth=1.0
        )
    ax.tick_params(labelsize=9, length=3, colors=ink)
    ax.spines[:].set_visible(False)
    fig.text(
        0.067, 0.465, "A full circle around z = +1", size=16, weight="bold", color=ink
    )
    fig.text(
        0.067,
        0.432,
        "Moving down shrinks the circle. Dotted lines mark the stable tangent.",
        size=10.5,
        color="#526675",
    )
    pol = fig.add_axes([0.16, 0.125, 0.75, 0.277])
    pol.imshow(
        local,
        extent=(0, 360, 80, 2),
        origin="upper",
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=3,
    )
    for angle in (
        stats["stable_tangent_degrees"],
        180 + stats["stable_tangent_degrees"],
    ):
        pol.axvline(angle, color="white", linestyle=(0, (2, 3)), lw=0.9)
    pol.set(
        xticks=[0, 90, 180, 270, 360],
        xticklabels=["0°", "90°", "180°", "270°", "360°"],
        yticks=[2, 20, 40, 60, 80],
        yticklabels=[
            r"$10^{-2}$",
            r"$10^{-20}$",
            r"$10^{-40}$",
            r"$10^{-60}$",
            r"$10^{-80}$",
        ],
        xlabel="Direction from +1",
        ylabel="Distance from +1",
    )
    pol.tick_params(labelsize=9, length=3, colors=ink)
    pol.spines[:].set_visible(False)
    for ax_ in (ax, pol):
        ax_.xaxis.label.set_color(ink)
        ax_.yaxis.label.set_color(ink)
    fig.text(
        0.067,
        0.049,
        "Rotate by 180° and shrink the radius by a factor of 1.22 × 10¹⁹.",
        size=10.5,
        color=ink,
    )
    total = stats["scaled_check"]["samples"]
    matched = total - stats["scaled_check"]["different_labels"]
    fig.text(
        0.067,
        0.023,
        f"{matched:,} of {total:,} tested pairs kept their basin color.",
        size=10.5,
        color=ink,
    )
    fig.savefig(PACKET / "basin-map.png", dpi=165, facecolor=paper)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["compute", "render"])
    parser.add_argument(
        "cache", type=Path, help="external .npz for the complete sampled grids"
    )
    args = parser.parse_args()
    if args.mode == "compute":
        if args.cache.resolve().is_relative_to(ROOT.resolve()):
            parser.error("Place the pixel cache outside the repository")
        compute(args.cache)
    else:
        render(args.cache)

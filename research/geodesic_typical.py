"""Compare parallel Jacobi fields over one physical lap of a seven-sheeted geodesic surface."""

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from flint import arb, ctx
from numba import njit
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
RAY = ROOT / "results/attractor_followup_2/geodesics/bounds.json"
N = 7
C = 3 / 8
LENGTH = 2 * np.pi * np.sqrt(99)


def area_bound(cells=4096):
    """Common projected area for 9/8 <= b < 11/8, rounded outward."""
    ctx.prec = 100
    h = arb.pi() / (2 * cells)
    total = arb(0)
    for j in range(cells):
        # cos(y) is positive and decreasing in this quadrant.
        lo = arb(11) / 8 * (j * h).cos() - arb(3) / 8
        hi = arb(9) / 8 * ((j + 1) * h).cos() + arb(3) / 8
        # Every cos(x) in [lo,hi] gives |cos(x)-b*cos(y)| <= c over this cell.
        if lo >= 1 or hi <= -1:
            continue
        lo = max(arb(-1), lo.upper())
        hi = min(arb(1), hi.lower())
        if lo < hi:
            total += h * (lo.acos() - hi.acos())
    # Four quadrants of y, two preimages of cos(x), fourteen sheets above (x,y).
    lower = (4 * 2 * 2 * N * total).lower()
    upper = (
        2 * arb.pi() * arb(99).sqrt() * (56 * arb.pi() / lower).sqrt() / arb(10).log()
    ).upper()
    assert lower > 116 and upper < arb(335) / 10
    return str(lower), str(upper)


def draw_starts(rng, b, count):
    """Uniform surface area by coordinate-line intersections, then uniform tangent angle."""
    starts = []
    while len(starts) < count:
        q = rng.uniform(0, 2 * np.pi, 3)
        axis = int(rng.integers(3))
        slot = int(rng.integers(2 * N))
        x, y, z = q
        if axis < 2:
            if slot >= 2:
                continue
            wave = (
                1 + C * np.exp(1j * (N * z - y))
                if axis == 0
                else b + C * np.exp(-1j * (x + N * z))
            )
            offset = b * np.cos(y) if axis == 0 else np.cos(x)
            value = -offset / abs(wave)
            if abs(value) >= 1:
                continue
            q[axis] = (-1) ** slot * np.arccos(value) - np.angle(wave)
        else:
            value = -(np.cos(x) + b * np.cos(y)) / C
            if abs(value) >= 1:
                continue
            q[2] = (
                (-1) ** slot * np.arccos(value) + 2 * np.pi * (slot // 2) - x + y
            ) / N
        q %= 2 * np.pi
        x, y, z = q
        sw = C * np.sin(x - y + N * z)
        g = np.array([-np.sin(x) - sw, -b * np.sin(y) + sw, -N * sw])
        normal = g / np.linalg.norm(g)
        # Coarea weights the proposal by |Nx|+|Ny|+|Nz| relative to surface area.
        if rng.uniform() > 1 / np.abs(normal).sum():
            continue
        e = np.cross(normal, np.eye(3)[np.argmin(abs(normal))])
        e /= np.linalg.norm(e)
        theta = rng.uniform(0, 2 * np.pi)
        v = e * np.cos(theta) + np.cross(normal, e) * np.sin(theta)
        starts.append(np.r_[q, v, 0.0, 0.0])
    return starts


@njit
def dynamics(t, state, b):
    x, y, z, vx, vy, vz, u, _logj = state
    w = x - y + N * z
    A = np.cos(x)
    B = b * np.cos(y)
    D = C * np.cos(w)
    sw = C * np.sin(w)
    gx = -np.sin(x) - sw
    gy = -b * np.sin(y) + sw
    gz = -N * sw
    gg = gx * gx + gy * gy + gz * gz
    normal_speed = gx * vx + gy * vy + gz * vz
    tx = vx - gx * normal_speed / gg
    ty = vy - gy * normal_speed / gg
    tz = vz - gz * normal_speed / gg
    vv = tx * tx + ty * ty + tz * tz
    hessian = -A * vx * vx - B * vy * vy - D * (vx - vy + N * vz) ** 2
    force = -(hessian + 3 * normal_speed + 2 * (A + B + D)) / gg
    k = (1 + b * b - C * C) / 2
    delta = (((1 + b) ** 2 - C * C) * (C * C - (b - 1) ** 2)) / 4
    minus_curvature = N * N * ((B + k * A) ** 2 + delta * A * A) / (gg * gg)
    return np.array(
        (
            vx,
            vy,
            vz,
            force * gx + (1 - vv) * tx,
            force * gy + (1 - vv) * ty,
            force * gz + (1 - vv) * tz,
            minus_curvature - u * u,
            u,
        )
    )


def integrate(job):
    power, initial, control = job
    b = 11 / 8 - 2.0 ** (-power)
    sol = solve_ivp(
        dynamics,
        (0, LENGTH),
        initial,
        args=(b,),
        method="DOP853",
        rtol=2e-12 if control else 2e-10,
        atol=1e-14 if control else 1e-12,
        max_step=0.06 if control else 0.14,
        t_eval=np.linspace(0, LENGTH, 257),
    )
    assert sol.success
    x, y, z, vx, vy, vz, u, logj = sol.y
    w = x - y + N * z
    level = np.max(abs(np.cos(x) + b * np.cos(y) + C * np.cos(w)))
    speed = np.max(abs(vx * vx + vy * vy + vz * vz - 1))
    tangent = np.max(
        abs(
            (-np.sin(x) - C * np.sin(w)) * vx
            + (-b * np.sin(y) + C * np.sin(w)) * vy
            - N * C * np.sin(w) * vz
        )
    )
    assert max(level, speed, tangent) < 1e-6 and u.min() > -1e-10
    return {
        "log10_parallel_growth": logj[-1] / np.log(10),
        "constraint": max(level, speed, tangent),
    }


def short_control(initial, b):
    """Direct second-order Jacobi integration with the ambient Hessian formula."""
    covectors = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, -1.0, N]])
    weights = np.array([1.0, b, C])

    def rhs(t, s):
        phase = covectors @ s[:3]
        g = -(weights * np.sin(phase)) @ covectors
        hessian = -covectors.T @ np.diag(weights * np.cos(phase)) @ covectors
        adj = np.stack(
            [
                np.cross(hessian[1], hessian[2]),
                np.cross(hessian[2], hessian[0]),
                np.cross(hessian[0], hessian[1]),
            ],
            axis=1,
        )
        gg = g @ g
        acceleration = -(s[3:6] @ hessian @ s[3:6]) / gg * g
        K = g @ adj @ g / gg**2
        return np.r_[s[3:6], acceleration, s[7], -K * s[6]]

    base = solve_ivp(
        dynamics,
        (0, 2),
        initial,
        args=(b,),
        method="DOP853",
        rtol=2e-11,
        atol=1e-13,
        max_step=0.04,
    )
    ref = solve_ivp(
        rhs,
        (0, 2),
        np.r_[initial[:6], 1.0, 0.0],
        method="Radau",
        rtol=2e-10,
        atol=1e-12,
        max_step=0.04,
    )
    assert base.success and ref.success
    discrepancy = float(abs(base.y[-1, -1] - np.log(ref.y[-2, -1])))
    assert discrepancy < 1e-7
    return discrepancy


def render(folder):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    packet = json.loads((folder / "growth.json").read_text())
    paper, ink, muted = "#fcfbf7", "#173440", "#52646c"
    teal, red = "#12837d", "#b94b28"
    fig = plt.figure(figsize=(9.1, 5.05), facecolor=paper)
    fig.text(
        0.071,
        0.921,
        "One lap, two very different scales of growth",
        fontsize=21,
        weight="bold",
        color=ink,
    )
    fig.text(
        0.073,
        0.844,
        "Initially parallel geodesics on the same surface, after 62.52 units of travel.",
        fontsize=12,
        color=muted,
    )
    ax = fig.add_axes([0.128, 0.19, 0.82, 0.55], facecolor=paper)
    levels = [r["power"] for r in packet["runs"]]
    limit = float(arb(packet["uniform_mean_log10_upper"]))
    ax.axvline(limit, c="#697a82", ls=(0, (4, 3)), lw=1.2, zorder=0)
    ax.text(
        limit + 1,
        4.38,
        "ceiling for mean log growth",
        fontsize=9,
        color=muted,
        ha="left",
    )
    rays = []
    jitter = np.random.default_rng(35)
    for j, r in enumerate(packet["runs"]):
        y = len(levels) - 1 - j
        obs = np.array(r["standard"] + r["tighter_independent"])
        ax.scatter(
            obs,
            y + jitter.uniform(-0.15, 0.15, len(obs)),
            s=11,
            facecolor=teal,
            edgecolors="none",
            alpha=0.55,
            zorder=2,
        )
        a, b = [float(arb(r[k])) for k in ("ray_log10_lower", "ray_log10_upper")]
        rays.append((a + b) / 2)
        ax.hlines(y, a, b, color=red, linewidth=3, zorder=4)
        ax.scatter([(a + b) / 2], [y], color=red, s=45, marker="D", zorder=5)
    ax.plot(rays, np.arange(len(levels))[::-1], color=red, alpha=0.5, lw=1.35, zorder=1)
    ax.set(
        yticks=np.arange(len(levels))[::-1],
        yticklabels=[f"$2^{{-{m}}}$" for m in levels],
        xlim=(0, 96),
        ylim=(-0.42, 4.62),
        xticks=np.arange(0, 91, 10),
    )
    ax.set_xlabel(
        "log$_{10}$ (final infinitesimal separation / initial separation)",
        labelpad=7,
        color=ink,
    )
    ax.set_ylabel("Distance to critical\nwave amplitude", labelpad=9, color=ink)
    ax.tick_params(axis="both", labelsize=10, colors=muted, length=3)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#96a5a5")
    fig.legend(
        handles=[
            Line2D(
                [],
                [],
                marker="o",
                ls="",
                markerfacecolor=teal,
                markeredgecolor="none",
                label="64 uniform launches per row",
            ),
            Line2D(
                [],
                [],
                marker="D",
                ls="",
                markerfacecolor=red,
                markeredgecolor=red,
                label="One closed straight ray",
            ),
        ],
        loc="upper left",
        bbox_to_anchor=(0.115, 0.078),
        ncol=2,
        frameon=False,
        fontsize=10,
        labelcolor=ink,
    )
    fig.savefig(folder / "parallel-rays.png", dpi=170, facecolor=paper)
    plt.close(fig)


def main(args):
    if args.render_only:
        return render(args.out)
    bounds = json.loads(RAY.read_text())
    area, upper = area_bound()
    rng = np.random.default_rng(84127)
    runs = []
    jobs = []
    for row in bounds:
        p = row["distance_power"]
        starts = draw_starts(rng, 11 / 8 - 2.0**-p, args.starts)
        # An independent ensemble guards against divergent pointwise comparisons in chaos.
        checks = draw_starts(rng, 11 / 8 - 2.0**-p, args.starts)
        for control, cohort in [(False, starts), (True, checks)]:
            jobs.extend((p, v, control) for v in cohort)
        # At a reflection-symmetric starting phase the transfer diagonal is half its trace.
        mu_low = arb(10) ** (arb(row["lower_hundredths"]) / 100)
        mu_hi = arb(10) ** (arb(row["upper_hundredths"]) / 100)
        J_low = ((mu_low + 1 / mu_low) / 2).log() / arb(10).log()
        J_hi = ((mu_hi + 1 / mu_hi) / 2).log() / arb(10).log()
        runs.append(
            {
                "power": p,
                "ray_log10_lower": str(J_low.lower()),
                "ray_log10_upper": str(J_hi.upper()),
                "short_direct_jacobi_error_in_log": short_control(
                    starts[0], 11 / 8 - 2.0 ** (-p)
                ),
            }
        )
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        values = list(pool.map(integrate, jobs, chunksize=1))
    for i, run in enumerate(runs):
        chunk = values[2 * i * args.starts : 2 * (i + 1) * args.starts]
        run["standard"] = [x["log10_parallel_growth"] for x in chunk[: args.starts]]
        run["tighter_independent"] = [
            x["log10_parallel_growth"] for x in chunk[args.starts :]
        ]
        run["max_sampled_constraint_error"] = max(x["constraint"] for x in chunk)
    args.out.mkdir(parents=True, exist_ok=True)
    result = {
        "seed": 84127,
        "divisions": 4096,
        "starts_per_cohort": args.starts,
        "lap_length": LENGTH,
        "area_lower": area,
        "uniform_mean_log10_upper": upper,
        "runs": runs,
    }
    (args.out / "growth.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    render(args.out)
    for run in runs:
        print(
            run["power"],
            np.mean(run["standard"]),
            np.mean(run["tighter_independent"]),
            flush=True,
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--starts", type=int, default=32)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--render-only", action="store_true", help="Render the saved growth.json."
    )
    main(p.parse_args())

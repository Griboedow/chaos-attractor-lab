"""Map three two-dimensional starting slices of the symmetric Sinkhorn example."""
import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from numba import njit, prange

OUT = Path(__file__).resolve().parents[1] / "results/sinkhorn-route"
COLORS = ["#f4f3eb", "#152b40", "#45c8ad", "#f69b76"]


@njit
def block(a, b):
    u, v = np.exp(a), np.exp(b)
    total = 1 + u + v
    reference = total + 7
    return np.log(reference/(total + 7*u)), np.log(reference/(total + 7*v))


@njit
def step(a, b, u, v):
    f, g = block(u, v)
    a, b = -.8*a + 1.8*f, -.8*b + 1.8*g
    f, g = block(a, b)
    return a, b, -.8*u + 1.8*f, -.8*v + 1.8*g


@njit(parallel=True)
def classify(grid, targets, slice_id, steps, tolerance):
    result = np.zeros((len(grid), len(grid)), dtype=np.uint8)
    for i in prange(len(grid)):
        for j in range(len(grid)):
            # Centered 3-log vector in two orthonormal coordinates; then anchor at row 1.
            X, Y = grid[j], grid[i]
            a = -np.sqrt(2.)*X
            b = -X/np.sqrt(2.) - np.sqrt(1.5)*Y
            if slice_id == 0:
                u, v = b-a, -a
            elif slice_id == 1:
                u, v = 0., 0.
            else:
                u, v = -b, a-b
            for t in range(steps):
                a, b, u, v = step(a, b, u, v)
                if (t+1) % 8 == 0:
                    for k in range(len(targets)):
                        delta = max(abs(a-targets[k,0]), abs(b-targets[k,1]),
                                    abs(u-targets[k,2]), abs(v-targets[k,3]))
                        if delta < tolerance:
                            result[i,j] = 1 if k == 0 else (2 if k <= 3 else 3)
                            break
                    if result[i,j]:
                        break
    return result


def render(grids, extent):
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.unicode_minus": True})
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 5.3), facecolor="#f7f7f2")
    fig.subplots_adjust(left=.048, right=.989, bottom=.205, top=.81, wspace=.115)
    fig.text(.043, .938, "One balanced matrix, three attractors", fontsize=24, color="#182b37")
    fig.text(.043, .873, "Diagonal entries 8, all others 1    ·    Target row and column sums 10    ·    Relaxation 9/5",
             fontsize=12, color="#385064")
    titles = ["Column start  (u₂, u₃, u₁)", "Column start  (1, 1, 1)", "Column start  (u₃, u₁, u₂)"]
    for i, (ax, img, title) in enumerate(zip(axes, grids, titles)):
        ax.imshow(img, origin="lower", extent=extent, interpolation="nearest",
                  cmap=ListedColormap(COLORS), vmin=-.5, vmax=3.5)
        ax.set_title(title, fontsize=12, color="#253d4d", loc="left", pad=7)
        ax.set_xticks([-30, 0, 30])
        ax.set_yticks([-30, 0, 30] if i == 0 else [])
        ax.tick_params(labelsize=9, colors="#546875", length=3)
        ax.spines[:].set_visible(False)
        if i == 0:
            ax.set_ylabel("Initial row-log contrast Y", fontsize=10, color="#445b69")
    fig.text(.515, .118, "Initial row-log contrast X", ha="center", fontsize=10, color="#445b69")
    handles = [Patch(facecolor=c, label=s) for c, s in zip(
        COLORS[1:], ["Balanced fixed point", "Largest row cycles 1 → 2 → 3", "Largest row cycles 1 → 3 → 2"])]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               labelcolor="#253d4d", bbox_to_anchor=(.5, .012), fontsize=10.5,
               handlelength=1.15, columnspacing=2.7)
    fig.savefig(OUT/"basins.png", dpi=145, facecolor=fig.get_facecolor())
    plt.close(fig)


def compare_direct_update():
    # Independent probability-domain formula, with no simplification for the chosen matrix.
    rng = np.random.default_rng(238919)
    kernel = 1 + 7*np.eye(3)
    error = 0.
    for state in rng.uniform(-60, 60, (200, 4)):
        u, v = np.exp(np.r_[0., state[:2]]), np.exp(np.r_[0., state[2:]])
        u = u**(-.8) * (10/(kernel @ v))**1.8
        v = v**(-.8) * (10/(kernel.T @ u))**1.8
        direct = np.r_[np.log(u[1:]/u[0]), np.log(v[1:]/v[0])]
        error = max(error, np.max(np.abs(direct - step(*state))))
    if error > 1e-12:
        raise ArithmeticError("independent matrix update disagrees")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="also repeat the common grid points at a tighter tolerance")
    args = parser.parse_args()
    phases = np.array(json.loads((OUT/"symmetric.json").read_text())["phase_seeds"])
    targets = np.vstack([np.zeros((1,4)), phases, phases[:,[1,0,3,2]]])
    grid = np.linspace(-32, 32, 513)
    tolerance, steps = 1e-8, 400
    maps = [classify(grid, targets, i, steps, tolerance) for i in range(3)]
    # Swapping the first two row labels reflects X and swaps both cycle handednesses.
    opposite = np.array([0,1,3,2])
    if not (np.array_equal(maps[0], opposite[maps[2][:,::-1]]) and
            np.array_equal(maps[1], opposite[maps[1][:,::-1]])):
        raise ArithmeticError("sampled basins violate reflection symmetry")
    if any((m == 0).any() for m in maps):
        raise ArithmeticError("unclassified grid points at the chosen horizon")
    if args.check:
        compare_direct_update()
        for i, m in enumerate(maps):
            tighter = classify(grid[::4], targets, i, 800, 1e-11)
            if not np.array_equal(tighter, m[::4, ::4]):
                raise ArithmeticError("classification changed on tighter check")
    half_pixel = (grid[1]-grid[0])/2
    render(maps, [grid[0]-half_pixel, grid[-1]+half_pixel]*2)
    summary = dict(interval=[-32,32], samples_per_axis=len(grid), max_steps=steps,
                   detection="Every eighth full update; sup-norm distance below 1e-8 to any target phase.",
                   columns=["unclassified","solution","row order 1-2-3","row order 1-3-2"],
                   counts=[np.bincount(m.ravel(), minlength=4).tolist() for m in maps])
    (OUT/"grid.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(summary)


if __name__ == "__main__":
    main()

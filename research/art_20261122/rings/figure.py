"""Render critical images and compare them to independent orbit ensembles."""
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/art_20261122/rings"
CASE = json.loads((ROOT / "results/wow_search_20261122/symmetry/annuli.json").read_text())["13745"]
_, lam, alpha, beta, gamma, omega, delta, e, h = map(float.fromhex, CASE["binary64_hex"])
L, C, E = complex(lam, omega), complex(alpha, delta), complex(e, h)
terms = [(L, 1, 0), (C, 3, 0), (beta/2, 18, 1),
         (beta/2, 18, -1), (gamma, 16, -1), (E, 33, -2)]
SOURCE = complex(*[float.fromhex(x) for x in json.loads((OUT / "cusp-witness.json").read_text())["source_binary64"]])


def step(z):
    return (L+C*abs(z)**2+beta*(z**17).real)*z + gamma*z.conjugate()**16 + E*z.conjugate()**33


def determinant(r, u):
    # F(r exp(i theta)) = exp(i theta) W(r,u), with u=17 theta.
    pieces = [c*r**d*np.exp(1j*k*u) for c, d, k in terms]
    radial = sum(v*d/r for v, (_, d, _) in zip(pieces, terms))
    angular = sum(1j*(1+17*k)*v for v, (_, _, k) in zip(pieces, terms))
    return (radial.conjugate()*angular).imag


def roots_at(angles, lo, hi):
    left = np.full(angles.shape, lo)
    right = np.full(angles.shape, hi)
    sign = np.sign(determinant(left, angles))
    assert np.all(sign*determinant(right, angles) < 0)
    for _ in range(49):
        mid = (left+right)/2
        same = sign*determinant(mid, angles) > 0
        left, right = np.where(same, mid, left), np.where(same, right, mid)
    return ((left+right)/2)*np.exp(1j*angles/17)


def branches():
    outer_edge = float(CASE["annuli"][1][1])
    start = brentq(lambda u: determinant(outer_edge, u), -.65, -.1)
    end = brentq(lambda u: determinant(outer_edge, u), -.1, .4)
    return (roots_at(np.linspace(0, 2*np.pi, 8193), .507, .509),
            roots_at(np.linspace(start, end, 4097), .775, .79))


TIP = step(SOURCE)
local_axis = step(SOURCE+1e-7)-step(SOURCE-1e-7)
AXIS = -local_axis/abs(local_axis)
X = np.array([-.0083, .0016])
Y = np.array([-.00048, .00048])


def population(seed, outer):
    # Independent starts, uniform in angle and radius in the inner trapping annulus.
    rng = np.random.default_rng(seed)
    z = rng.uniform(.36, .62, 2048)*np.exp(2j*np.pi*rng.random(2048))
    counts = np.zeros((100, 210))
    strip_counts = np.zeros((2, 2), dtype=int)
    folded_arc = (step(outer)-TIP)/AXIS
    turn = np.argmax(folded_arc.real)
    sides = (folded_arc[:turn+1], folded_arc[turn:][::-1])
    sides = sorted(sides, key=lambda s: np.interp(-.005, s.real, s.imag))
    for _ in range(1024):
        z = step(z)
    for _ in range(16):
        batch = []
        for _ in range(256):
            z = step(z)
            batch.append(z)
        sample = np.concatenate(batch)
        u = ((np.angle(sample)-np.angle(TIP)+np.pi/17) % (2*np.pi/17))-np.pi/17
        local = (abs(sample)*np.exp(1j*(u+np.angle(TIP)))-TIP)/AXIS
        counts += np.histogram2d(local.imag, local.real, bins=counts.shape, range=[Y, X])[0]
        active = local[(local.real > -.0075) & (local.real < -.002)]
        for i, side in enumerate(sides):
            # Positive normal here points into the narrow tongue on both sides.
            distance = (active.imag - np.interp(active.real, side.real, side.imag))*(1-2*i)
            strip_counts[i] += [np.count_nonzero((distance > 0) & (distance < .00002)),
                                np.count_nonzero((distance < 0) & (distance > -.00002))]
    return counts, strip_counts


def make_figure(inner, outer, counts):
    bg = "#101d28"
    pale = "#e9dfcf"
    gold = "#f4b36f"
    fig = plt.figure(figsize=(12.5, 7), facecolor=bg)
    ax = fig.add_axes([.012, .045, .54, .92], facecolor=bg, aspect="equal")
    palette = LinearSegmentedColormap.from_list("age", ["#fff0b0", "#f1a675", "#a96786", "#71a9b0", "#6798b4"])
    for kind, curve in enumerate((inner, outer)):
        stages = []
        for k in range(26):
            curve = step(curve)
            stages.append(curve)
        for k in reversed(range(len(stages))):
            color = palette(k/(len(stages)-1))
            for turn in range(17):
                line = stages[k]*np.exp(2j*np.pi*turn/17)
                ax.plot(line.real, line.imag, color=color, alpha=.67 if k < 10 else .35,
                        lw=(.45 if kind == 0 else .62), solid_capstyle="round")
    ax.plot(TIP.real, TIP.imag, marker="o", ms=7, mfc="none", mec="#f7f6e7",
            mew=1, zorder=5)
    ax.set(xlim=(-.82, .82), ylim=(-.82, .82))
    ax.axis("off")
    fig.text(.565, .905, "A cusp inside the crown", fontsize=24, color=pale, fontfamily="serif")
    fig.text(.567, .853, "Two critical curves, traced through 26 map updates", color="#a7bac1", fontsize=10)
    density = fig.add_axes([.575, .455, .388, .29], facecolor=bg)
    cmap = LinearSegmentedColormap.from_list("density", [bg, "#274652", "#58858a", "#a1c5bb", "#f2d9ac"])
    density.imshow(np.log1p(counts), extent=[*X, *Y], origin="lower", cmap=cmap, aspect="auto",
                   interpolation="nearest")
    density.set_title("Visits near one cusp", color=pale, loc="left", fontsize=12, pad=8)
    theory = fig.add_axes([.575, .105, .388, .25], facecolor=bg, sharex=density, sharey=density)
    profile = (step(outer)-TIP)/AXIS
    theory.plot(profile.real, profile.imag, color=gold, lw=1.8)
    theory.set_title("Image of the critical arc", color=pale, loc="left", fontsize=12, pad=8)
    for panel in (density, theory):
        panel.set(xlim=X, ylim=Y, yticks=[-.0004, 0, .0004])
        panel.tick_params(colors="#9cb3bd", labelsize=8, length=3)
        panel.spines[["top", "right"]].set_visible(False)
        panel.spines[["left", "bottom"]].set_color("#586b76")
        panel.ticklabel_format(style="sci", axis="y", scilimits=(-4, -4), useMathText=True)
        panel.yaxis.get_offset_text().set(color="#9cb3bd", fontsize=8)
    density.tick_params(labelbottom=False)
    theory.set_xlabel("Distance along the cusp tangent", fontsize=9, color="#b5c3c8", labelpad=5)
    fig.savefig(OUT / "postcritical-crown.png", dpi=158, facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    inner, outer = branches()
    counts, first = population(13745, outer)
    _, second = population(5109, outer)
    make_figure(inner, outer, counts)
    print("Visits displayed:", int(counts.sum()))
    for seed, count in ((13745, first), (5109, second)):
        print(f"Seed {seed}, inside/outside counts (lower, upper): {count.tolist()}.")

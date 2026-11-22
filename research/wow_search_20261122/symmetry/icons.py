"""Scout and inspect cyclically equivariant polynomial maps from physical orbit points."""

from pathlib import Path
import argparse
import math
import os

import numpy as np
from numba import njit, prange, set_num_threads
from PIL import Image, ImageDraw
from scipy.stats import qmc


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/wow_search_20261122/symmetry"
CACHE = Path(
    os.environ.get(
        "WOW_SYMMETRY_CACHE",
        "D:/Codex/chaos-attractor-cache/wow-search-20261122/symmetry",
    )
)
CACHE.mkdir(parents=True, exist_ok=True)
SEED = 2026112207
INITIAL = complex(0.247, 0.193)


@njit
def power(z, n):
    t = 1.0 + 0.0j
    while n:
        if n & 1:
            t *= z
        z *= z
        n //= 2
    return t


@njit
def step(z, p):
    n = int(p[0])
    q = power(z, n - 2)
    zn1 = q * z
    zn = zn1 * z
    u = z.real * z.real + z.imag * z.imag
    a = complex(p[2], p[6])
    c = complex(p[7], p[8])
    b = p[3]
    g = p[4]
    h = complex(p[1], p[5]) + a * u + b * zn.real
    f = h * z + g * zn1.conjugate() + c * (zn1 * zn).conjugate()
    # Wirtinger derivatives with respect to z and its conjugate.
    A = h + a * u + (0.5 * b * n) * zn
    B = (
        a * z * z
        + (0.5 * b * n) * z * zn1.conjugate()
        + g * (n - 1) * q.conjugate()
        + c * (2 * n - 1) * (zn1 * zn1).conjugate()
    )
    return f, A, B


@njit
def orbit(p, initial, burn, count, windows):
    z = initial
    v = 1.0 + 0.0j
    seq = np.empty(count, np.complex128)
    rates = np.zeros((windows, 2))
    for i in range(burn + count):
        f, A, B = step(z, p)
        if (
            abs(f) > 4.0
            or abs(f) < 1.0e-15
            or not np.isfinite(f.real)
            or not np.isfinite(f.imag)
        ):
            return seq, rates, False
        w = A * v + B * v.conjugate()
        s = abs(w)
        if s < 1.0e-100 or not np.isfinite(s):
            return seq, rates, False
        v = w / s
        z = f
        if i >= burn:
            k = i - burn
            seq[k] = z
            j = k * windows // count
            rates[j, 0] += math.log(s)
            det = A.real * A.real + A.imag * A.imag - B.real * B.real - B.imag * B.imag
            rates[j, 1] += math.log(max(abs(det), 1.0e-300)) - math.log(s)
    return seq, rates * (windows / count), True


def parameters(power_of_two):
    u = qmc.Sobol(10, scramble=True, seed=SEED).random_base2(power_of_two)
    anchors = [
        (5, -1.806, 1.806, 0, 1, 0),
        (6, -2.7, 5.0, 2.0, 1.0, 0),
        (16, 2.39, -2.5, -0.1, 0.9, -0.15),
        (3, 1.5, -1.0, -0.1, -0.8, 0),
        (9, -2.4, 3.8, 0.85, 1.0, 0),
        (4, -2.35, 3.1, 0.2, 1.0, 0.02),
    ]
    rows = np.empty((len(u), 9))
    for i, r in enumerate(u):
        n, l, a, b, g, w = anchors[min(int(r[0] * len(anchors)), len(anchors) - 1)]
        n = max(3, min(20, n + int(np.floor(7 * r[1])) - 3))
        rows[i] = (
            n,
            l + (r[2] - 0.5) * 1.24,
            a + (r[3] - 0.5) * 1.75,
            b + (r[4] - 0.5) * 1.7,
            g + (r[5] - 0.5) * 1.1,
            w + (r[6] - 0.5) * 0.63,
            (r[7] - 0.5) * 1.25,
            (r[8] - 0.5) * 1.3,
            (r[9] - 0.5) * 1.3,
        )
    return rows


@njit(parallel=True)
def batch(rows):
    result = np.zeros((len(rows), 7))
    for i in prange(len(rows)):
        p = rows[i]
        z, r, ok = orbit(p, INITIAL, 1200, 3072, 4)
        if not ok or r[:, 0].min() < 0.006:
            continue
        rad = np.abs(z)
        R = rad.max()
        if R < 0.01:
            continue
        sectors = np.zeros(18)
        bins = np.zeros((48, 48))
        ang = np.exp(1j * int(p[0]) * np.angle(z))
        shape = np.abs(np.mean(ang))
        for j in range(len(z)):
            t = z[j] / R
            k = int((math.atan2(t.imag, t.real) + math.pi) / (2 * math.pi) * 18) % 18
            sectors[k] += 1
            x = min(47, max(0, int((t.real + 1) * 24)))
            y = min(47, max(0, int((t.imag + 1) * 24)))
            bins[y, x] += 1
        # Broad angular coverage plus sparse but resolved traces.
        occupied = np.sum(bins > 0)
        if sectors.min() < 8 or occupied < 85:
            continue
        result[i] = np.array(
            (
                r[:, 0].mean(),
                r[:, 1].mean(),
                r[:, 0].min(),
                occupied,
                sectors.min() / len(z),
                shape,
                np.min(rad) / R,
            )
        )
    return result


def score_gallery(rows, metrics, maximum):
    idx = np.flatnonzero(metrics[:, 0])
    print("Scout flags", len(idx), flush=True)
    candidates = []
    for i in idx:
        m = metrics[i]
        rank = (
            m[5]
            + 0.20 * min(1.0, m[0] / 0.19)
            + 0.40 * (1.0 - min(1.0, m[3] / 1450.0))
            + 0.25 * min(1.0, m[4] * 40)
        )
        candidates.append((rank, i))
    candidates.sort(reverse=True)
    shortlist = []
    # Mix sparse filaments, radial contrast, and strongly anisotropic sectors.
    orders = [
        candidates,
        sorted(candidates, key=lambda x: metrics[x[1], 5], reverse=True),
        sorted(
            candidates,
            key=lambda x: abs(np.log(metrics[x[1], 3] / 420)) - 0.22 * metrics[x[1], 6],
        ),
        sorted(
            candidates,
            key=lambda x: abs(np.log(metrics[x[1], 3] / 950)) - 0.35 * metrics[x[1], 5],
        ),
    ]
    for k in range(len(idx)):
        for order in orders:
            i = order[k][1]
            if i not in shortlist:
                shortlist.append(i)
        if len(shortlist) >= maximum:
            break
    for j in range(0, len(shortlist), 40):
        tile = Image.new("RGB", (8 * 165, 5 * 184), (7, 10, 19))
        draw = ImageDraw.Draw(tile)
        for l, i in enumerate(shortlist[j : j + 40]):
            seq, r, ok = orbit(rows[i], INITIAL, 12000, 65536, 8)
            if not ok:
                continue
            im = portrait(seq, 160)
            a, b = (l % 8) * 165, (l // 8) * 184
            tile.paste(im, (a, b))
            draw.text((a + 4, b + 160), f"{i}  n={int(rows[i, 0])}", fill="white")
        tile.save(CACHE / f"gallery-{j // 40:02d}.jpg", quality=92)
    return shortlist


def hist(z, size):
    radius = np.abs(z).max() * 1.035
    d = np.histogram2d(
        z.imag, z.real, bins=size, range=[(-radius, radius), (-radius, radius)]
    )[0][::-1]
    return d, radius


def portrait(z, size):
    d, _ = hist(z, size)
    t = np.log1p(d) / np.log1p(max(2, np.quantile(d[d > 0], 0.997)))
    knots = [0, 0.12, 0.28, 0.46, 0.69, 0.85, 1.0]
    colors = np.array(
        [
            [5, 9, 19],
            [20, 44, 77],
            [35, 96, 134],
            [42, 158, 167],
            [141, 205, 162],
            [244, 191, 113],
            [255, 237, 200],
        ]
    )
    rgb = np.stack(
        [np.interp(np.minimum(t, 1), knots, colors[:, i]) for i in range(3)], axis=-1
    ).astype("uint8")
    return Image.fromarray(rgb)


def parameter_for(case, power_of_two=14):
    if case >= 200000:
        from odd import parameters_odd, BASE

        return parameters_odd()[case - BASE]
    return parameters(power_of_two)[case]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["scout", "render"])
    parser.add_argument("--power", type=int, default=14)
    parser.add_argument("--count", type=int, default=2000000)
    parser.add_argument("--case", type=int, default=0)
    parser.add_argument("--gallery", type=int, default=160)
    parser.add_argument("--threads", type=int, default=3)
    args = parser.parse_args()
    set_num_threads(args.threads)
    rows = parameters(args.power)
    if args.action == "scout":
        metrics = batch(rows)
        shortlist = score_gallery(rows, metrics, args.gallery)
        np.savez_compressed(
            CACHE / f"scout-{args.power}.npz",
            parameters=rows,
            metrics=metrics,
            shortlist=shortlist,
        )
        print("shortlist", shortlist, flush=True)
    else:
        p = parameter_for(args.case, args.power)
        seq, r, ok = orbit(p, INITIAL, 100000, args.count, 20)
        assert ok
        name = f"cyclic-{args.case}"
        np.savez_compressed(
            CACHE / f"{name}.npz",
            points=np.column_stack((seq.real, seq.imag)),
            rates=r,
            p=p,
        )
        portrait(seq, 1200).save(OUT / f"{name}.png")
        print(name, "p", p.tolist(), "LEs", r.mean(axis=0).tolist(), flush=True)


if __name__ == "__main__":
    main()

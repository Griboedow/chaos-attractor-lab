"""Compare two ordered torus orbits through one sevenfold observation."""
import argparse
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/art_20261122/rosette"
MOD = 2**61 - 1
START = (15485863, 32452843)
N = 2**24
BLOCK = 2048
RADIUS = 2
WIDTH = 1024
WEIGHT = 0.27
SHIFT = (sqrt(2)/100, sqrt(3)/100)


def matrix_power(n):
    """Exact modular power of the cat map, using Python integer products."""
    a = np.array([[2, 1], [1, 1]], dtype=object)
    product = np.eye(2, dtype=object)
    while n:
        if n & 1:
            product = (product @ a) % MOD
        a = (a @ a) % MOD
        n //= 2
    return product


def phase_blocks(system):
    # Each column is a disjoint consecutive stretch of the one archived orbit.
    indices = np.arange(N//BLOCK) * BLOCK
    starts = []
    state = np.array(START, dtype=object)
    if system == "cat":
        jump = matrix_power(BLOCK)
        for _ in indices:
            starts.append(state)
            state = (jump @ state) % MOD
        x, y = np.array(starts, dtype=np.int64).T
    for step in range(BLOCK+1):
        if system == "cat":
            yield x.astype(float)/MOD, y.astype(float)/MOD
            x, y = (2*x+y) % MOD, (x+y) % MOD
        else:
            t = indices+step
            yield (START[0]/MOD+t*SHIFT[0]) % 1, (START[1]/MOD+t*SHIFT[1]) % 1


def observation(x, y):
    u = np.exp(2j*pi*x)
    return u**3 + (2/3)*u.conj()**4 + WEIGHT*np.exp(2j*pi*y)


def pixels(z):
    x = ((z.real+RADIUS)*WIDTH/(2*RADIUS)).astype(int)
    y = ((RADIUS-z.imag)*WIDTH/(2*RADIUS)).astype(int)
    return x + WIDTH*y


def measure(system, collect_pairs):
    full = np.zeros(WIDTH**2, dtype=np.uint32)
    halfway = np.zeros_like(full)
    phases = iter(phase_blocks(system))
    z = observation(*next(phases))
    cross = 0j
    pairs = []
    for x, y in phases:
        w = observation(x, y)
        if collect_pairs:
            # The first columns are consecutive blocks of the original orbit.
            ncols = 2**20//BLOCK
            pairs.append(np.stack((z[:ncols], w[:ncols])))
        np.add.at(full, pixels(z), 1)
        np.add.at(halfway, pixels((z+w)/2), 1)
        cross += np.vdot(z, w)
        z = w
    assert int(full.sum()) == int(halfway.sum()) == N
    pair_array = np.concatenate(pairs, axis=1) if pairs else None
    return full.reshape(WIDTH, WIDTH), halfway.reshape(WIDTH, WIDTH), cross/N, pair_array


def coarse(hist, width):
    factor = WIDTH//width
    return hist.reshape(width, factor, width, factor).sum((1, 3)).astype(np.uint32)


def image(hist):
    from PIL import Image
    # The gradient depends only on the measured number of visits per pixel.
    palette = np.array([[12, 12, 29], [40, 28, 68], [84, 46, 111],
                        [162, 72, 124], [226, 130, 119], [255, 220, 158]])
    anchors = [0, .15, .36, .57, .79, 1.]
    v = np.minimum(1, np.arcsinh(hist/12)/np.arcsinh(24.))
    rgb = np.stack([np.interp(v, anchors, channel) for channel in palette.T], axis=-1)
    return Image.fromarray(rgb.astype(np.uint8))


def film(sources):
    from PIL import Image, ImageDraw, ImageFont
    from matplotlib import get_data_path
    width = 352
    font = ImageFont.truetype(str(Path(get_data_path())/"fonts/ttf/DejaVuSans.ttf"), 17)
    frames = []
    for t in np.linspace(0, 1, 41):
        frame = Image.new("RGB", (2*width, width+36), (12, 12, 29))
        draw = ImageDraw.Draw(frame)
        for side, ((z, w), label) in enumerate(zip(sources, ("Rigid rotation", "Cat map"))):
            points = (1-t)*z + t*w
            x = ((points.real+RADIUS)*width/(2*RADIUS)).astype(int)
            y = ((RADIUS-points.imag)*width/(2*RADIUS)).astype(int)
            counts = np.bincount(x+width*y, minlength=width**2).reshape(width, width)
            # Match the reference density per unit area in the still.
            light = counts * (N*width**2/(len(z)*WIDTH**2))
            frame.paste(image(light), (side*width, 32))
            draw.text(((side+.5)*width, 15), label, fill="#f3d7be", font=font, anchor="mm")
        frames.append(frame)
    durations = [110]*len(frames)
    durations[0], durations[20] = 950, 650
    frames[0].save(OUT/"between-samples.webp", save_all=True, append_images=frames[1:],
                   duration=durations, loop=0, quality=82, method=3)


def cusps():
    from flint import arb, ctx
    from fractions import Fraction
    ctx.prec = 144
    # Isolate a*kappa(c)=1 by exact rational bisection of its squared equation.
    lo, hi = Fraction(1, 2), Fraction(1)
    for _ in range(100):
        c = (lo+hi)/2
        lhs = (81*(72*c-13))**2
        rhs = 100**2*(145-144*c)**3
        if lhs > rhs:
            hi = c
        else:
            lo = c
    def ball(f):
        return arb(f.numerator)/f.denominator
    angle = ball(lo).union(ball(hi)).acos()
    corners = []
    for j in range(7):
        for sign in (1, -1):
            u = (2*j*arb.pi()+sign*angle)/7
            x = (3*u).cos()+2*(4*u).cos()/3
            y = (3*u).sin()-2*(4*u).sin()/3
            dx = -3*(3*u).sin()-8*(4*u).sin()/3
            dy = 3*(3*u).cos()-8*(4*u).cos()/3
            scale = (arb(27)/100)/(dx**2+dy**2).sqrt()
            corners.append((x-scale*dy, y+scale*dx))
    gaps = []
    for i, (x, y) in enumerate(corners):
        for u, v in corners[:i]:
            dx, dy = x-u, y-v
            gaps.append((dx*dx + dy*dy).sqrt())
    assert all(gap > arb("0.47") for gap in gaps)
    return min(float(gap.lower()) for gap in gaps)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--render", action="store_true")
    parser.add_argument("--film", action="store_true", help="Render chord positions for consecutive samples")
    args = parser.parse_args()
    cat, cat_mid, corr_cat, cat_pairs = measure("cat", args.film)
    rot, rot_mid, corr_rot, rot_pairs = measure("rotation", args.film)
    var = 1 + (2/3)**2 + WEIGHT**2
    alpha, beta = 2*pi*np.array(SHIFT)
    corr_exact = np.exp(3j*alpha) + (2/3)**2*np.exp(-4j*alpha) + WEIGHT**2*np.exp(1j*beta)
    assert abs(corr_cat)/var < .004 and abs(corr_rot-corr_exact)/var < .0002
    tv = np.abs(coarse(cat, 64).astype(float)-coarse(rot, 64)).sum()/(2*N)
    assert tv < .025
    receipt = dict(samples=N, bin_grid=WIDTH, histogram_tv_64=float(tv),
                   cat_lag_one=float(abs(corr_cat)/var), rotation_lag_one=float(abs(corr_rot)/var),
                   rotation_exact_lag_one=float(abs(corr_exact)/var), cusp_separation=cusps())
    print(json.dumps(receipt, indent=2))
    if args.render:
        import matplotlib.pyplot as plt
        OUT.mkdir(parents=True, exist_ok=True)
        image(cat).save(OUT/"sevenfold-rosette.png", optimize=True)
        fig, axes = plt.subplots(1, 2, figsize=(9, 4.8), facecolor="#0c0c1d", layout="constrained")
        for ax, hist, name in zip(axes, (rot_mid, cat_mid), ("Rigid rotation", "Cat map")):
            ax.imshow(image(hist), interpolation="bilinear")
            ax.set_title(name, color="#f3d7be", fontsize=13)
            ax.axis("off")
        fig.savefig(OUT/"midpoints.png", dpi=134, facecolor=fig.get_facecolor())
        plt.close(fig)
        for name, hist in (("cat", cat), ("rotation", rot)):
            np.savez_compressed(OUT/(name+"-counts.npz"), counts=coarse(hist, 256))
        (OUT/"measurements.json").write_text(json.dumps(receipt,indent=2)+"\n", encoding="utf8",newline="\n")
    else:
        saved = json.loads((OUT/"measurements.json").read_text())
        for name in receipt:
            assert np.isclose(receipt[name], saved[name], rtol=1e-7, atol=1e-10), name
        for name, hist in (("cat", cat), ("rotation", rot)):
            with np.load(OUT/(name+"-counts.npz")) as data:
                assert np.array_equal(data["counts"], coarse(hist, 256))
    if args.film:
        film((rot_pairs, cat_pairs))


if __name__ == "__main__":
    main()

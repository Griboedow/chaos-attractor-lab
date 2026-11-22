"""Render Newton basin colors on their nine-scale logarithmic quotient torus."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numba import get_num_threads, njit, prange, set_num_threads
from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/art_20261122/newton"
EXACT = ROOT / "results/wow_search_20261122/boundary/certificate.json"
TAU = 2*math.pi
Q = 5/6
PAL = ["639cda", "73c8b5", "445b9e", "e2b078", "46a9bf", "b087b2",
       "ffca72", "ff9556", "ed5847", "dc477a", "b556a7", "8a6ad6",
       "4f9bc3", "69c9ad", "ddce70"]
COLORS = np.array([list(bytes.fromhex(s)) for s in PAL], dtype=float)/255


def inputs():
    data = json.loads(EXACT.read_text())
    def value(row):
        re, im = row["center_dyadic"]
        return complex(math.ldexp(int(re[0]), re[1]), math.ldexp(int(im[0]), im[1]))
    return value(data["parameter_certificate"]), np.array([value(r) for r in data["root_certificates"]])


@njit
def texture(angle, scale, a, roots):
    cell = int(math.floor(scale))
    fraction = scale - cell
    w = 1000 * Q**(-fraction) * np.exp(-1j*angle)
    for step in range(350):
        if abs(w) < 2:
            for idx in range(6):
                dist = abs(w-roots[idx])
                if dist < 1e-7:
                    fine = math.log(math.log(max(dist, 1e-16))/math.log(1e-7))/math.log(2)
                    return idx, step-fine-fraction
            if abs(w) < 1e-7:
                return 6+(-step-cell-1)%9, step-fraction
        if abs(w) > 3:
            w = Q*w + a*(Q*w-1)/(6*w**5-a)
        else:
            w = (5*w**6-a)/(6*w**5-a)
    return -1, 350.


@njit(parallel=True)
def render(width, height, a, roots):
    # Orthographic view of a circular torus. Bounds and camera do not affect the coloring.
    rgb = np.empty((height, width, 3), dtype=np.uint8)
    lab = np.full((height, width), -2, dtype=np.int16)
    incl = 0.61
    si, co = math.sin(incl), math.cos(incl)
    R, r = 1.43, .61
    for py in prange(height):
        y = (.5-(py+.5)/height)*3.52
        for px in range(width):
            x = ((px+.5)-width/2)/height*3.52
            yy, zz = y*si-4.5*co, y*co+4.5*si
            bg = max(0., 1.-((x/3.4)**2+((y+.1)/2.3)**2))
            col = np.array([.029+.033*bg, .051+.047*bg, .078+.060*bg])
            t = 0.
            for _ in range(170):
                Y, Z = yy+t*co, zz-t*si
                rho = math.sqrt(x*x+Y*Y)
                dist = math.hypot(rho-R, Z)-r
                if dist < .000018 or t > 8.:
                    break
                t += dist
            if t <= 8. and dist < .000018:
                u = math.atan2(Y, x)
                v = math.atan2(Z, rho-R)
                # Logarithmic scale wraps around the large circle. Physical angle wraps around the tube.
                scale = (9*(u/TAU + .24)) % 9.
                azimuth = .5*(v+u+.45)
                angle = (-math.pi-.01 + 2*math.atan2(.055*math.sin(azimuth), math.cos(azimuth))) % TAU
                label, time = texture(angle, scale, a, roots)
                lab[py, px] = label
                cv, sv = math.cos(v), math.sin(v)
                nx, ny, nz = math.cos(u)*cv, math.sin(u)*cv, sv
                light = max(0., -.48*nx-.45*ny+.75*nz)
                back = max(0., .57*nx+.56*ny+.60*nz)
                front = max(0., -co*ny+si*nz)
                spec = max(0., -.26*nx-.76*ny+.595*nz)**38
                edge = (1.-front)**3
                if label < 0:
                    base = np.array([.1,.13,.19])
                else:
                    base = COLORS[label].copy()
                # Slow capture darkens the filaments. The corrected capture time repeats every scale.
                detail = max(.28, min(1.05, math.exp(-(time-44)/42)))
                etched = .86 + .14*(.5 + .5*math.cos(TAU*time))**6
                depth = .34+.74*light+.18*back+.17*front
                col = base*depth*detail*etched + np.array([.30,.39,.44])*(.38*edge+.75*spec)
                if label >= 6:
                    col += base*.12
            for c in range(3):
                rgb[py, px, c] = int(255*min(1.,max(0.,col[c])))
    return rgb, lab


@njit(parallel=True)
def strip(tile_width, height, a, roots):
    # Repeat a single cell by the exact label action. Every row holds the linearizing angle fixed.
    rgb = np.empty((height, 9*tile_width, 3), dtype=np.uint8)
    phase_counts = np.zeros((height, 9), dtype=np.int64)
    for iy in prange(height):
        angle = -math.pi-.01 + ((iy+.5)/height-.5)*.085
        for ix in range(tile_width):
            s = (ix+.5)/tile_width
            label, time = texture(angle, s, a, roots)
            detail = max(.35, min(1.06, math.exp(-(time-44)/58)))
            etched = .90+.10*(.5+.5*math.cos(TAU*time))**6
            for k in range(9):
                current = label if label < 6 else 6+(label-6-k)%9
                color = COLORS[current] if current >= 0 else np.array([.08,.10,.14])
                if current >= 6:
                    phase_counts[iy, current-6] += 1
                for channel in range(3):
                    rgb[iy, k*tile_width+ix, channel] = min(255, int(255*color[channel]*detail*etched))
    return rgb, phase_counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--width", type=int, default=2100)
    ap.add_argument("--height", type=int, default=1360)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--supersample", type=int, default=2)
    args = ap.parse_args()
    set_num_threads(min(args.workers, get_num_threads()))
    a, roots = inputs()
    scale = args.supersample
    pixels, labels = render(args.width*scale, args.height*scale, a, roots)
    image = Image.fromarray(pixels)
    if scale > 1:
        image = image.resize((args.width, args.height), Image.Resampling.LANCZOS)
    image.save(OUT / "nine-scale-newton-torus.png", optimize=True)
    active = labels != -2
    counts = np.bincount(labels[active]+1, minlength=16).tolist()
    report = {"case": 6902, "size": [args.width, args.height], "supersample_per_axis": scale,
              "surface_subpixels": int(active.sum()), "unresolved": counts[0],
              "outcome_subpixels_in_label_order": counts[1:], "palette": PAL,
              "coordinates": "scale=9*(u/(2*pi)+0.24) mod 9; angle=-pi-.01+2*atan2(.055*sin(b),cos(b)) mod 2pi; b=(v+u+.45)/2"}
    strip_rgb, phase_counts = strip(260, 560, a, roots)
    assert np.all(phase_counts == phase_counts[:, [0]])
    Image.fromarray(strip_rgb).save(OUT / "nine-scale-weave.png", optimize=True)
    report["weave"] = {"shape": [560, 2340], "angle_span_radians": .085,
                       "angle_center": "-pi-.01", "per_phase_pixels": phase_counts.sum(axis=0).tolist(),
                       "rows_with_equal_phase_counts": 560}
    (OUT/"render.json").write_text(json.dumps(report, indent=2)+"\n")
    print(report)


if __name__ == "__main__":
    main()

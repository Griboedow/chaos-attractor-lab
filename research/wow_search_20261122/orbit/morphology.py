"""Inspect data-derived silhouettes with a deterministic diversity selection."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import label
from PIL import Image, ImageDraw


def coordinates(z):
    origin = z.mean(0)
    s, v = np.linalg.eigh(np.cov(z.T))
    axes = v[:, ::-1]
    for j in range(3):
        if axes[np.argmax(abs(axes[:, j])), j] < 0:
            axes[:, j] *= -1
    if np.linalg.det(axes) < 0:
        axes[:, 2] *= -1
    return (z - origin) @ axes, axes


def feature(z):
    zz, axes = coordinates(z)
    rms = np.sqrt((zz * zz).sum(1).mean())
    zz /= max(rms, 1e-8)
    sig = zz.std(0)
    features = [*(sig[1:] / max(sig[0], 1e-12))]
    for norm in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0.8, 0.6]):
        n = np.asarray(norm, dtype=float)
        a = np.cross(n, [1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0])
        a /= np.sqrt(a @ a)
        b = np.cross(n, a)
        uv = np.stack((zz @ a, zz @ b), 1)
        q = np.quantile(uv, [0.005, 0.995], axis=0)
        box = np.maximum(q[1] - q[0], 1e-7)
        ij = np.floor((uv - q[0]) / box * 47).astype(int).clip(0, 47)
        freq = np.bincount(ij[:, 0] * 48 + ij[:, 1], minlength=2304).reshape(48, 48)
        features.extend(
            [np.mean(freq > 0), np.mean(freq >= 4), label(freq >= 3)[1] / 60]
        )
    step = zz[1:] - zz[:-1]
    features.extend(
        np.quantile(np.linalg.norm(step, axis=1), [0.05, 0.5, 0.9]).tolist()
    )
    lagstats = np.stack(
        [
            np.quantile(np.linalg.norm(zz[k:] - zz[:-k], axis=1), [0.10, 0.50])
            for k in range(1, 25)
        ]
    )
    features.extend(
        [
            np.min(lagstats[:, 0]),
            np.min(lagstats[:, 1]),
            np.argmin(lagstats[:, 0]) / 24,
            np.argmin(lagstats[:, 1]) / 24,
        ]
    )
    zs = z[:, 2]
    width = np.quantile(zs, 0.98) - np.quantile(zs, 0.02)
    f = np.histogram(zs, 24)[0]
    features.extend(
        [
            np.sum((f[1:-1] > f[:-2]) & (f[1:-1] > f[2:])) / 10,
            (zs.max() - zs.min()) / max(width, 1e-6),
        ]
    )
    return features


def project(z, mode=0):
    zz, _ = coordinates(z)
    if mode == 1:
        return z[:, [0, 1]], z[:, 2]
    if mode == 2:
        return zz[:, [0, 2]], zz[:, 1]
    if mode == 3:
        return z[:, [0, 2]], z[:, 1]
    x = zz[:, 0] * 0.91 + zz[:, 2] * 0.415
    y = zz[:, 1] * 0.94 - zz[:, 2] * 0.34
    return np.column_stack((x, y)), zz[:, 2]


def tile(z, size=(156, 130), mode=0):
    w, h = size
    xy, color = project(z, mode)
    mid = (xy.max(0) + xy.min(0)) / 2
    span = np.ptp(xy, axis=0)
    scale = 0.94 * min((w - 2) / max(span[0], 1e-9), (h - 2) / max(span[1], 1e-9))
    ix = np.minimum(
        w - 1,
        np.maximum(0, np.rint((xy[:, 0] - mid[0]) * scale + (w - 1) / 2).astype(int)),
    )
    iy = np.minimum(
        h - 1,
        np.maximum(0, np.rint(-(xy[:, 1] - mid[1]) * scale + (h - 1) / 2).astype(int)),
    )
    idx = iy * w + ix
    counts = np.bincount(idx, minlength=w * h).reshape(h, w)
    u = (color - np.quantile(color, 0.01)) / max(
        np.quantile(color, 0.99) - np.quantile(color, 0.01), 1e-9
    )
    cs = np.bincount(idx, weights=np.clip(u, 0, 1), minlength=w * h).reshape(
        h, w
    ) / np.maximum(counts, 1)
    strength = np.minimum(1, np.log1p(counts) / np.log(4))
    col = np.array([4, 12, 25]) + strength[:, :, None] * (
        (1 - cs[:, :, None]) * np.array([64, 206, 210])
        + cs[:, :, None] * np.array([253, 161, 83])
    )
    return Image.fromarray(np.asarray(col, dtype="uint8"))


def main(a):
    paths = sorted(a.cache.glob("deep-*.npz"))
    metrics = []
    f = []
    ids = []
    params = []
    clouds = []
    for path in paths:
        if path.name == "deep-design.npz":
            continue
        with np.load(path) as q:
            mm = q["metrics"]
            zz = q["points"]
            pp = q["params"]
            ii = q["ids"]
            for i in np.flatnonzero(mm[:, :4].min(1) > 0.01):
                z = zz[i].astype(float)
                ff = feature(z)
                if not np.isfinite(ff).all():
                    continue
                clouds.append(zz[i])
                f.append(ff)
                metrics.append(mm[i])
                ids.append(int(ii[i]))
                params.append(pp[i])
    f = np.asarray(f)
    mm = np.asarray(metrics)
    np.savez_compressed(
        a.cache / "morphology.npz", ids=ids, params=params, metrics=mm, features=f
    )
    lo, hi = np.quantile(f, [0.08, 0.92], axis=0)
    v = np.clip((f - lo) / np.maximum(hi - lo, 0.03), -2, 4)
    seeds = np.unique(
        np.concatenate([np.argmax(v, 0), np.argmin(v, 0), [np.argmax(mm[:, 14])]])
    ).tolist()
    dist = np.min([np.linalg.norm(v - v[k], axis=1) for k in seeds], axis=0)
    while len(seeds) < min(144, len(f)):
        dist[seeds] = -1
        i = np.argmax(dist)
        seeds.append(int(i))
        dist = np.minimum(dist, np.linalg.norm(v - v[i], axis=1))
    for page in range((len(seeds) + 35) // 36):
        canvas = Image.new("RGB", (6 * 162, 6 * 156), (25, 33, 44))
        draw = ImageDraw.Draw(canvas)
        for j, ix in enumerate(seeds[page * 36 : (page + 1) * 36]):
            row, col = divmod(j, 6)
            canvas.paste(tile(clouds[ix]), (col * 162 + 3, row * 156 + 2))
            draw.text(
                (col * 162 + 5, row * 156 + 133),
                f"{ids[ix]}    {mm[ix, :4].mean():.3f}",
                fill="white",
            )
        canvas.save(a.cache / f"morphology-{page}.jpg", quality=93)
    (a.cache / "visual-ids.json").write_text(json.dumps([ids[i] for i in seeds]))
    print(len(ids), len(seeds), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--cache", type=Path, required=True)
    main(p.parse_args())

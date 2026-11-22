"""Broad low-discrepancy screen and native-coordinate morphology shortlist."""

import argparse
import json
from pathlib import Path
import numpy as np
import pyopencl as cl
from scipy.stats import qmc


def designs(power):
    u = qmc.Sobol(6, scramble=True, seed=229866505).random_base2(power)
    return np.ascontiguousarray(
        np.column_stack(
            (
                0.25 + 2.5 * u[:, 0],
                0.44 + 0.545 * u[:, 1],
                (2 * u[:, 2] - 1) * np.pi,
                (2 * u[:, 3] - 1) * 19,
                (2 * u[:, 4] - 1) * 27,
                0.0025 * (0.85 / 0.0025) ** u[:, 5],
            )
        ),
        dtype="f4",
    )


def gpu():
    dev = next(
        d
        for p in cl.get_platforms()
        for d in p.get_devices()
        if d.type & cl.device_type.GPU
    )
    ctx = cl.Context([dev])
    return (
        ctx,
        cl.CommandQueue(ctx),
        cl.Program(ctx, Path(__file__).with_name("survey.cl").read_text())
        .build()
        .survey,
    )


def run(env, p, burn, win, keep):
    ctx, q, k = env
    p = np.ascontiguousarray(p, "f4")
    n = len(p)
    x = np.empty((n, 16), "f4")
    z = np.empty((max(1, keep), n, 4), "f4")
    flags = cl.mem_flags
    pb = cl.Buffer(ctx, flags.READ_ONLY | flags.COPY_HOST_PTR, hostbuf=p)
    xb = cl.Buffer(ctx, flags.WRITE_ONLY, x.nbytes)
    zb = cl.Buffer(ctx, flags.WRITE_ONLY, z.nbytes)
    k(
        q,
        (n,),
        None,
        pb,
        xb,
        zb,
        np.int32(n),
        np.int32(burn),
        np.int32(win),
        np.int32(keep),
    )
    cl.enqueue_copy(q, x, xb).wait()
    if keep:
        cl.enqueue_copy(q, z, zb).wait()
    return x, z[:, :, :3]


def features(p, m):
    n = len(m)
    cov = np.zeros((n, 3, 3))
    for i in range(3):
        cov[:, i, i] = m[:, 7 + i]
    for j, (a, b) in enumerate([(0, 1), (0, 2), (1, 2)]):
        cov[:, a, b] = cov[:, b, a] = m[:, 10 + j]
    ev = np.maximum(np.linalg.eigvalsh(cov), 1e-15)
    ev /= ev.sum(1)[:, None]
    f = np.column_stack(
        (
            m[:, :4].mean(1),
            m[:, :4].std(1),
            ev[:, :2],
            np.log1p(m[:, 13]),
            m[:, 14],
            m[:, 15],
            np.log(p[:, 5]),
            np.sqrt(np.maximum(m[:, 9], 0.0))
            / np.sqrt(np.maximum(m[:, 7] + m[:, 8], 1e-7)),
        )
    )
    return ev, np.nan_to_num(f)


def main(a):
    a.cache.mkdir(parents=True, exist_ok=True)
    env = gpu()
    p = designs(a.power)
    fx = a.cache / "screen.npz"
    if fx.exists():
        m = np.load(fx)["metrics"]
    else:
        m = np.concatenate(
            [run(env, p[i : i + 4096], 3584, 512, 0)[0] for i in range(0, len(p), 4096)]
        )
        np.savez_compressed(fx, params=p, metrics=m)
    ev, f = features(p, m)
    good = np.flatnonzero((m[:, :4].min(1) > 0.012) & (ev[:, 0] > 0.0003))
    nt = min(2400, len(good))
    z = f[good]
    width = np.maximum(
        np.quantile(z, 0.90, axis=0) - np.quantile(z, 0.10, axis=0), 0.001
    )
    z = z / width
    sel = [
        np.argmax(np.min(m[good, :4], axis=1)),
        np.argmax(ev[good, 0]),
        np.argmax(m[good, 14]),
    ]
    dist = np.min([np.linalg.norm(z - z[j], axis=1) for j in sel], axis=0)
    for _ in range(nt - len(sel)):
        dist[sel] = -1
        j = np.argmax(dist)
        sel.append(j)
        dist = np.minimum(dist, np.linalg.norm(z - z[j], axis=1))
    ids = good[sel]
    np.savez_compressed(a.cache / "deep-design.npz", ids=ids, params=p[ids])
    for i in range(0, len(ids), 128):
        out = a.cache / f"deep-{i:04}.npz"
        if out.exists():
            continue
        small = ids[i : i + 128]
        mm, cc = run(env, p[small], 24576, 4096, 4096)
        np.savez_compressed(
            out, ids=small, params=p[small], metrics=mm, points=cc.transpose(1, 0, 2)
        )
    print(
        json.dumps(
            dict(
                unique=len(p),
                short_positive=int(np.sum(m[:, :4].min(1) > 0.012)),
                positive_with_full_geometry=len(good),
                deep=len(ids),
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True, type=Path)
    ap.add_argument("--power", type=int, default=17)
    main(ap.parse_args())

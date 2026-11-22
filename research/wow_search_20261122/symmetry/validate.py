"""Finite-time diagnostics at fixed independent starts, with an independent polar evaluator."""

import hashlib
import json
import math

import numpy as np
from numba import njit
from icons import CACHE, OUT, INITIAL, parameters, parameter_for, orbit, step

CASES = (1723, 13745, 14725, 11450, 5529, 205188, 202080)
STARTS = (
    INITIAL,
    -0.31 + 0.11j,
    0.433 - 0.173j,
    0.64 + 0.07j,
    -0.216 + 0.36j,
    -0.72 - 0.08j,
    0.93 + 0.04j,
)
COUNT = 2000000
SHORT_COUNT = 400000
BURN = 100000
WINDOWS = 20


@njit
def polar(z, p):
    n = int(p[0])
    r = abs(z)
    theta = math.atan2(z.imag, z.real)
    phi = n * theta
    c, s = math.cos(phi), math.sin(phi)
    sn, cn = math.sin(2 * phi), math.cos(2 * phi)
    re = (
        r * (p[1] + p[2] * r * r + p[3] * r**n * c)
        + p[4] * r ** (n - 1) * c
        + r ** (2 * n - 1) * (p[7] * cn + p[8] * sn)
    )
    im = (
        r * (p[5] + p[6] * r * r)
        - p[4] * r ** (n - 1) * s
        + r ** (2 * n - 1) * (p[8] * cn - p[7] * sn)
    )
    return complex(re, im) * complex(math.cos(theta), math.sin(theta))


@njit
def replay(p, z, burn, count):
    data = np.empty(count, np.complex128)
    for k in range(burn + count):
        z = polar(z, p)
        if abs(z) > 4.0 or not np.isfinite(z.real):
            return data, False
        if k >= burn:
            data[k - burn] = z
    return data, True


def normalize(x):
    return x / x.sum()


def histogram(z, extent):
    counts = np.histogram2d(
        z.real, z.imag, bins=64, range=[(-extent, extent), (-extent, extent)]
    )[0]
    assert int(counts.sum()) == len(z)
    return normalize(counts)


def sectors(z, n):
    return normalize(
        np.histogram(np.angle(z), bins=n, range=(-np.pi, np.pi))[0]
    ).tolist()


def summary(x):
    return dict(
        mean=x.mean(0).tolist(), minimum=x.min(0).tolist(), maximum=x.max(0).tolist()
    )


def periodic_test(z):
    # Maximum late sup-norm error for each candidate lag, over 1,024 consecutive states.
    block = z[-1024:]
    values = [
        np.maximum(
            abs((block - z[-1024 - k : -k]).real), abs((block - z[-1024 - k : -k]).imag)
        ).max()
        for k in range(1, 513)
    ]
    return dict(
        minimum_over_lags=float(min(values)), best_lag=int(np.argmin(values)) + 1
    )


def main():
    scout = np.load(CACHE / "scout-14.npz")
    assert np.array_equal(scout["parameters"], parameters(14))
    out = {
        "seed": 2026112207,
        "scout_unique": len(scout["parameters"]),
        "scout_flags": int(np.count_nonzero(scout["metrics"][:, 0])),
        "gallery_examined": len(scout["shortlist"]),
        "selected": list(CASES),
        "parameter_order": [
            "n",
            "lambda",
            "alpha",
            "beta",
            "gamma",
            "omega",
            "delta",
            "Re(eta)",
            "Im(eta)",
        ],
        "start_order": [[v.real, v.imag] for v in STARTS],
        "warmup": BURN,
        "portrait_observed": COUNT,
        "windows": WINDOWS,
        "alternate_observed": SHORT_COUNT,
        "cases": {},
    }
    odd = np.load(CACHE / "odd-scout.npz")
    all_rows = np.concatenate((scout["parameters"], odd["parameters"]))
    assert len(np.unique(all_rows, axis=0)) == len(all_rows)
    out["odd_scout_unique"] = len(odd["metrics"])
    out["odd_scout_flags"] = int(np.count_nonzero(odd["metrics"][:, 0]))
    for i in CASES:
        f = CACHE / f"cyclic-{i}.npz"
        with np.load(f) as packet:
            xy = packet["points"]
            z = xy[:, 0] + 1j * xy[:, 1]
            p = packet["p"]
            rate = packet["rates"]
        assert np.array_equal(p, parameter_for(i)) and len(z) == COUNT
        R = float(abs(z).max() * 1.04)
        main_h = histogram(z, R)
        half_tv = (
            0.5
            * abs(histogram(z[: COUNT // 2], R) - histogram(z[COUNT // 2 :], R)).sum()
        )
        checks = []
        for initial in STARTS[1:]:
            seq, qr, ok = orbit(p, initial, BURN, SHORT_COUNT, 20)
            row = dict(start=[initial.real, initial.imag], complete=bool(ok))
            if ok:
                if abs(seq).max() <= R:
                    row["histogram_tv_to_portrait"] = float(
                        0.5 * abs(histogram(seq, R) - main_h).sum()
                    )
                row.update(
                    rates=summary(qr),
                    radius_range=[float(abs(seq).min()), float(abs(seq).max())],
                    sector_fractions=sectors(seq, int(p[0])),
                )
            checks.append(row)
        polar_z, ok = replay(p, INITIAL, BURN, SHORT_COUNT)
        assert ok and abs(polar_z).max() <= R
        errors = []
        for t in z[:: len(z) // 64]:
            t = complex(t)
            u, A, B = step(t, p)
            direct = polar(t, p)
            fd = []
            for v in [1.0, 1j]:
                approx = (polar(t + v * 1.0e-6, p) - polar(t - v * 1.0e-6, p)) / (
                    2.0e-6
                )
                fd.append(abs(approx - A * v - B * complex(v).conjugate()))
            errors.append([abs(u - direct), *fd])
        image_file = OUT / f"cyclic-{i}.png"
        assert image_file.is_file()
        out["cases"][str(i)] = dict(
            parameters=p.tolist(),
            parameters_hex=[float(x).hex() for x in p],
            rates=summary(rate),
            sector_fractions=sectors(z, int(p[0])),
            radius_range=[float(abs(z).min()), float(abs(z).max())],
            tv_first_vs_second_half=float(half_tv),
            other_starts=checks,
            polar_replay_tv=float(0.5 * abs(histogram(polar_z, R) - main_h).sum()),
            largest_single_step_polar_discrepancy_and_fd_column_errors=np.max(
                errors, axis=0
            ).tolist(),
            lag_check=periodic_test(z),
            external=f.name,
            external_cache="wow-search-20261122/symmetry/" + f.name,
            cloud_sha256_float64_le=hashlib.sha256(
                xy.astype("<f8", copy=False).tobytes()
            ).hexdigest(),
            npz_sha256=hashlib.file_digest(f.open("rb"), "sha256").hexdigest(),
            png_sha256=hashlib.file_digest(image_file.open("rb"), "sha256").hexdigest(),
        )
        line = out["cases"][str(i)]
        print(
            i,
            "LE",
            line["rates"],
            "sector range",
            [min(line["sector_fractions"]), max(line["sector_fractions"])],
            "TV",
            round(half_tv, 4),
            round(line["polar_replay_tv"], 4),
            "alt",
            [
                (
                    v["complete"],
                    round(v.get("histogram_tv_to_portrait", -1), 4),
                    v.get("rates", {}).get("mean"),
                )
                for v in checks
            ],
            flush=True,
        )
    (OUT / "numerics.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()

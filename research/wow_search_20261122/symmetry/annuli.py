"""Propose radial trapping annuli by sampled polar bounds, then replay with Arb."""

import json

import numpy as np
from icons import CACHE, INITIAL, orbit, parameters, parameter_for


def envelope(p, lo, hi):
    n, l, a, b, g, w, d, e, h = p
    n = int(n)
    radius = np.linspace(lo, hi, 161)[:, None]
    phase = np.arange(513) * 2 * np.pi / 512
    t = np.exp(-1j * phase)
    v = (
        radius * (l + 1j * w + (a + 1j * d) * radius**2 + b * radius**n * t.real)
        + g * radius ** (n - 1) * t
        + (e + 1j * h) * radius ** (2 * n - 1) * t * t
    )
    m = abs(v)
    return m.min(), m.max()


def propose(i, p=None, q=None, pad=0.0001):
    if p is None:
        p = parameters(14)[i]
    seq, rate, ok = orbit(p, INITIAL, 100000, 120000, 20)
    assert ok
    r = abs(seq)
    if q is None:
        q = {5529: 4, 15960: 5}.get(i, 2)
    initial = np.array([(r[k::q].min(), r[k::q].max()) for k in range(q)])
    bands = initial.copy()

    for t in range(140):
        old = bands.copy()
        for k in range(q):
            lo, hi = envelope(p, *old[k])
            j = (k + 1) % q
            bands[j] = [min(bands[j, 0], lo - pad), max(bands[j, 1], hi + pad)]
        if min(bands[:, 0]) <= 0 or max(bands[:, 1]) > np.max(r) * 1.09:
            print(i, "escaped", t, bands, flush=True)
            return
        if np.max(abs(bands - old)) < 1e-11:
            break
    print(i, "candidate", t, np.round(initial, 6), np.round(bands, 6), flush=True)
    return np.round(bands, 7)


if __name__ == "__main__":
    for case, period, pad in [
        (5529, 4, 0.0001),
        (13745, 2, 0.0001),
        (205188, 3, 0.0001),
        (202080, 5, 0.000013),
    ]:
        bands = propose(case, parameter_for(case), period, pad)
        if bands is None:
            raise RuntimeError(f"No proposed bands for case {case}")
        ordered = bands[np.argsort(bands[:, 0])]
        assert np.all(ordered[1:, 0] > ordered[:-1, 1])
        (CACHE / f"proposed-bands-{case}.json").write_text(
            json.dumps(bands.tolist(), indent=2)
        )

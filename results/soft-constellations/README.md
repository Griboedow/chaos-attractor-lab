# Gaussian-posterior circle experiment

The [scientific report](../../docs/soft-constellations.md) describes the
temperature sweep, selected numerical regimes and limits of the evidence.
The [priority review](../../docs/soft-constellation-priority.md) identifies
known mean-shift, MMSE, proximal and softmax formulations.

Nine completed run directories are inventoried: scout, all-flag follow-up,
long runs, dense rendering, 33-start cohorts, two precision settings,
reference-atlas comparison and standalone-averaging/alternating-lobe controls.
`inventory.json` lists every one of their 241 files, including partial
reports and exact executed source snapshots. The 178 optional arrays total
664,774,925 bytes and remain in the hash-addressed D: cache and original
run directories. Exact inputs and compact metrics are retained under
distinct paths, as recorded in `omissions.json`; see the
[storage instructions](../../docs/data-storage.md).

The exact-model parameters, including every site's hexadecimal binary64
coordinates, are in `runs/soft-constellation-long01/report.json`. Screen
identifiers are local to this study. In particular, Case 156 here is the
positive-temperature continuation of old finite-point Case 1725, not old
finite-point Case 156.

## Reproduction

Use new external directories. The original finite-point survey supplies
the six input geometries unchanged.

```text
python research/search_raar_constellations.py --out D:/Codex/new-circle-scout --power 9
python research/search_soft_constellations.py --source D:/Codex/new-circle-scout/screen.npz --out D:/Codex/new-soft-scout
python research/follow_soft_constellations.py --source D:/Codex/new-soft-scout/screen.npz --out D:/Codex/new-soft-follow
python research/extend_soft_constellations.py --source D:/Codex/new-soft-scout/screen.npz --out D:/Codex/new-soft-long --ids 27 29 59 92 122 156 191 --selection "Recorded visual selection and two-expanding-direction comparison"
python research/render_soft_constellations.py --source D:/Codex/new-soft-follow --out D:/Codex/new-soft-render --images docs/assets/experimental/new-soft-replay --ids 29 59 92 122 156 191
python research/cohort_soft_constellations.py --source D:/Codex/new-soft-long --out D:/Codex/new-soft-cohort --ids 29 156
python research/check_soft_constellation_precision.py --source D:/Codex/new-soft-scout/screen.npz --out D:/Codex/new-soft-precision70 --ids 29 156 --digits 70
python research/check_soft_constellation_precision.py --source D:/Codex/new-soft-scout/screen.npz --out D:/Codex/new-soft-precision100 --ids 29 156 --digits 100
python research/compare_constellation_catalog.py --source D:/Codex/new-soft-long --out D:/Codex/new-soft-catalog --ids 29 156
python research/analyze_soft_constellations.py --source D:/Codex/new-soft-long --cohort D:/Codex/new-soft-cohort --dense D:/Codex/new-soft-render --out D:/Codex/new-soft-structure
```

The precision control requires mpmath 1.3.0. Its optional `--mpmath-path`
appends an external library directory without overriding NumPy/Numba.
The recorded environment uses Python 3.12.14, NumPy 2.5.4, Numba 0.68.1,
SciPy 1.18.1 and Matplotlib 3.11.2. Bitwise trajectory identity is an
environment-specific replay check, not a promise across platforms.

After restoring the archived optional inputs, the accounting/replay command is

```text
python research/verify_soft_campaign.py --out D:/Codex/new-soft-verification
```

This checks all original hashes, the six hard controls, all seven
five-million-step runs, all six two-million-state renderings and their
histograms. It separately compares the tangent calculation with LAPACK QR
and determinant accumulation over 2,048 steps per case. The [completed receipt](verification/report.json) records a successful replay.
These are numerical checks, not interval certificates.

The shared [portrait catalog](../../records/portrait-catalog.json) supplies
the paired main-page rendering and 24-second rotation. Each rotation panel
uses the last 250,000 consecutive states of its complete rendering input;
the data and physical scales stay fixed throughout. First drawing dates
are recorded in the [date registry](../../records/render-dates.json).

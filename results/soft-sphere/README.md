# Gaussian posterior mean and sphere reflection

The [numerical article](../../docs/soft-sphere/manuscript.md) reports the
complete selection protocol, spatial candidates, competing regular outcomes,
equivalent formulations and unresolved novelty. The main portraits are
Cases 155 and 251 **of this sphere-temperature sweep**.

Ten completed run directories contain 262 inventoried files. The 123 optional
arrays total 487,169,399 bytes and remain in the hash-addressed external
cache and original D: directories. Compact inputs and summaries have distinct
filenames recorded in `omissions.json`; full trajectories were not silently
replaced by thinned data. See the [storage instructions](../../docs/data-storage.md).

The inventory includes all 429 screen cases, all 138 flags and follow-ups,
all seven selected long runs, five dense renderings, four 33-start cohorts,
four cases at each of two precision settings, the 108-regime comparison,
standalone averaging and regular-outcome controls, and numerical verification.
Source snapshots are the bytes actually executed. Exact hexadecimal model
parameters and sites are in
[the verification report](runs/soft-sphere-verification01/report.json).

## Reproduction

Run into new external directories. The first two commands reproduce the
earlier hard-circle and hard-sphere input surveys without changing them.

```text
python research/search_raar_constellations.py --out D:/Codex/new-circle --power 9
python research/search_spherical_constellations.py --source D:/Codex/new-circle/screen.npz --out D:/Codex/new-sphere
python research/search_soft_sphere.py --mode screen --source D:/Codex/new-sphere/screen.npz --out D:/Codex/new-soft-sphere-scout
python research/search_soft_sphere.py --mode follow --source D:/Codex/new-soft-sphere-scout/screen.npz --out D:/Codex/new-soft-sphere-follow
python research/search_soft_sphere.py --mode follow --source D:/Codex/new-soft-sphere-scout/screen.npz --out D:/Codex/new-soft-sphere-long --ids 155 196 250 251 253 388 389 --count 5000000 --warm 100000 --selection "Recorded visual shortlist after reviewing all twelve pages of the complete follow-up"
python research/render_soft_sphere.py --source D:/Codex/new-soft-sphere-follow --out D:/Codex/new-soft-sphere-dense --images docs/assets/experimental/new-soft-sphere-replay --ids 155 196 251 388 389
python research/cohort_soft_sphere.py --source D:/Codex/new-soft-sphere-long --out D:/Codex/new-soft-sphere-cohort --ids 155 196 251 389
python research/check_soft_sphere_precision.py --source D:/Codex/new-soft-sphere-scout/screen.npz --out D:/Codex/new-soft-sphere-precision70 --ids 155 196 251 389 --digits 70
python research/check_soft_sphere_precision.py --source D:/Codex/new-soft-sphere-scout/screen.npz --out D:/Codex/new-soft-sphere-precision100 --ids 155 196 251 389 --digits 100
python research/compare_constellation_catalog.py --source D:/Codex/new-soft-sphere-long --out D:/Codex/new-soft-sphere-catalog --ids 155 196 251 389
python research/analyze_soft_sphere.py --source D:/Codex/new-soft-sphere-long --cohort D:/Codex/new-soft-sphere-cohort --dense D:/Codex/new-soft-sphere-dense --out D:/Codex/new-soft-sphere-structure
python research/verify_soft_sphere.py --scout D:/Codex/new-soft-sphere-scout --follow D:/Codex/new-soft-sphere-follow --long D:/Codex/new-soft-sphere-long --dense D:/Codex/new-soft-sphere-dense --out D:/Codex/new-soft-sphere-verification
```

The environment used Python 3.12.14, NumPy 2.5.4, Numba 0.68.1, SciPy
1.18.1, Matplotlib 3.11.2 and mpmath 1.3.0. The precision tool's optional
`--mpmath-path` appends an external library directory without replacing the
NumPy/Numba environment. Bitwise replay is environment-specific.

The successful verification repeats every screen output, accounts for all
follow-up flags and return tests, repeats the seven five-million-step
observations and all five two-million-state renderings, and reconstructs
their histograms. Its local comparisons use an independently expressed map,
finite differences, LAPACK QR and determinant accumulation. These are
numerical controls, not interval enclosures or a proof of chaotic attraction.

The [portrait catalog](../../records/portrait-catalog.json) and shared
renderers produce the main paired image and fixed-data 24-second rotation.
Each animated panel uses the last 250,000 consecutive rendering states.
The [drawing registry](../../records/render-dates.json) records actual
first-render dates. The twelve preview pages contain every flagged case;
their legacy `leading_xy` metadata key is the full 3D leading rate here.

The frozen reference-atlas identifier `ZhouChen` is a software label with
unverified bibliographic identification; see the article's source caveat.
Reference-cloud distances are not historical novelty certificates.

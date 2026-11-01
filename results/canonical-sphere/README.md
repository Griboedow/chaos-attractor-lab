# Canonical posterior-reflection campaign

The [numerical article](../../docs/canonical-sphere/manuscript.md) reports
the folded ribbon of Case 163 and spatial network of Case 653, along with
the complete selection history and limitations. Identifiers are local to
this campaign. Neither chaotic attraction nor historical uniqueness is
certified.

Fifteen completed run directories contain 417 inventoried files. The 271
optional arrays total 888,372,437 bytes and remain in the external cache
and original D: run directories. Compact exact inputs and small summaries
are retained under distinct filenames; full trajectories were not silently
thinned. See [data storage](../../docs/data-storage.md), `inventory.json`
and `omissions.json`.

The packet contains all 656 screen cases, all 262 follow-ups, eight long
observations, seven dense renderings, four 33-start leading-rate cohorts,
two 33-start full-spectrum cohorts, independent precision controls,
periodic-capture tests, a 108-regime geometric comparison, and numerical
replay. Executed source snapshots accompany each run. The
[strict JSON summary](selected-summary.json) provides the selected parameters
and diagnostics; exact hexadecimal parameters and sites are in the two
verification reports.

The original follow-up serializer retains unused legacy metadata
`base_id=-1` and `relative_tau=NaN`. The latter is a Python JSON extension,
not a valid strict-JSON number. Those bytes are preserved for provenance;
neither field enters the map. `selected-summary.json` omits these two
undefined metadata fields and is strict JSON. Physical temperatures are
the finite `tau` values in the complete scout protocol.

## Reproduction

Use new external output directories. Obtain the earlier hard-sphere and
soft-sphere screen inputs by following the
[preceding packet's reproduction instructions](../soft-sphere/README.md).
Their hashes are recorded in this campaign's scout protocol. The commands
below use those newly reproduced input paths.

```text
python research/canonical_soft_sphere.py --hard-source D:/Codex/new-sphere/screen.npz --soft-source D:/Codex/new-soft-sphere-scout/screen.npz --out D:/Codex/new-canonical-scout
python research/search_soft_sphere.py --mode follow --source D:/Codex/new-canonical-scout/screen.npz --out D:/Codex/new-canonical-follow
python research/render_projection_gallery.py --source D:/Codex/new-canonical-follow --out docs/assets/experimental/new-canonical-preview --title "Canonical posterior sphere: every follow-up"
python research/search_soft_sphere.py --mode follow --source D:/Codex/new-canonical-scout/screen.npz --out D:/Codex/new-canonical-long --ids 7 163 361 407 565 638 653 --warm 100000 --count 5000000 --selection "Recorded visual shortlist after reviewing all 22 preview pages"
python research/search_soft_sphere.py --mode follow --source D:/Codex/new-canonical-scout/screen.npz --out D:/Codex/new-canonical-simplelong --ids 148 --warm 100000 --count 5000000 --selection "Explicit added control: simplest temperature 0.06 near Case 7"
python research/render_soft_sphere.py --source D:/Codex/new-canonical-long --out D:/Codex/new-canonical-dense --images docs/assets/experimental/new-canonical-dense --ids 7 163 361 407 565 638 653 --count 2000000
python research/cohort_soft_sphere.py --source D:/Codex/new-canonical-long --out D:/Codex/new-canonical-cohort --ids 7 163 361 653
python research/check_soft_sphere_precision.py --source D:/Codex/new-canonical-scout/screen.npz --out D:/Codex/new-canonical-precision70 --ids 7 163 653 --digits 70
python research/check_soft_sphere_precision.py --source D:/Codex/new-canonical-scout/screen.npz --out D:/Codex/new-canonical-precision100 --ids 7 163 653 --digits 100
python research/check_sphere_full_precision.py --source D:/Codex/new-canonical-scout/screen.npz --out D:/Codex/new-canonical-fullprecision --ids 653 --digits 70 100 --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages
python research/check_sphere_full_precision.py --source D:/Codex/new-canonical-scout/screen.npz --out D:/Codex/new-canonical-simpleprecision --ids 148 --digits 70 100 --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages
python research/inspect_canonical_sphere.py --source D:/Codex/new-canonical-scout/screen.npz --follow D:/Codex/new-canonical-follow --out D:/Codex/new-canonical-inspection --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages
python research/analyze_canonical_sphere.py --scout D:/Codex/new-canonical-scout --cohort D:/Codex/new-canonical-cohort --out D:/Codex/new-canonical-structure
python research/compare_constellation_catalog.py --source D:/Codex/new-canonical-long --out D:/Codex/new-canonical-catalog --ids 163 653
python research/verify_soft_sphere.py --scout D:/Codex/new-canonical-scout --follow D:/Codex/new-canonical-follow --long D:/Codex/new-canonical-long --dense D:/Codex/new-canonical-dense --out D:/Codex/new-canonical-verification
python research/verify_soft_sphere.py --scout D:/Codex/new-canonical-scout --follow D:/Codex/new-canonical-follow --long D:/Codex/new-canonical-simplelong --dense D:/Codex/new-canonical-dense --out D:/Codex/new-canonical-simpleverification
```

The recorded environment is Python 3.12.14, NumPy 2.5.4, Numba 0.68.1,
SciPy 1.18.1, Matplotlib 3.11.2 and mpmath 1.3.0. Supply an installed
mpmath directory to the precision scripts' path argument. Bitwise replay
depends on the environment. The comparison additionally requires the 108
frozen reference trajectories identified by its report; use the data
registry to obtain their exact cached inputs. No reference equations are
silently substituted.

The successful replay checks every scout output and follow-up return test,
all eight five-million-step observations, and all seven two-million-state
renderings with reconstructed histograms. Independent local map expressions,
finite differences, LAPACK QR and accumulated log determinants provide
arithmetic controls. The 70/100-digit finite-separation spectra use no
production Jacobian. These are numerical tests, not interval proofs.

The [portrait catalog](../../records/portrait-catalog.json) and shared
renderers produce the main image and 24-second rotation. Each animated
panel uses the last 250,000 consecutive rendering states. All 22 preview
pages and all seven dense images were visually inspected. The preview's
legacy `leading_xy` field is the full three-dimensional leading rate here.
The [drawing registry](../../records/render-dates.json) records first renders.

The [literature review receipt](literature-review.json) distinguishes full
equation sections, abstracts, inspected figures, and search-only leads.
Its scope is bounded: an unsuccessful exact-formula search is not evidence
that no equivalent attractor has ever been published.

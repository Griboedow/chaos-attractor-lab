# Higher-order and nonuniform product-posterior searches

This frozen packet extends the [product-posterior article](../../docs/product-alphabet/manuscript.md)
without changing its established posterior/reflection operator. A regular
five-, seven- and nine-level survey gives mostly thin frames and ribbons.
A second survey with nonuniform three-, four- and five-level alphabets
gives the additional folded-shell numerical Cases 6 and 7.

**Identifiers are local to each survey.** In both surveys, Cases 0 and 1
reproduce the earlier product Cases 988 and 489. A higher-order Case 7 and
a nonuniform Case 7 would not denote the same parameters. Exact alphabet
entries are in each scout protocol and compact input records.

## Evidence and selection

| Stage | Regular higher-order survey | Nonuniform survey |
| --- | --- | --- |
| Declared inputs | 851: two controls, 81 grid cases, 256 Sobol geometries at three orders | 770: two controls, 256 Sobol geometries at three orders |
| Sobol design | Eight dimensions, seed 2026110152 | Seventeen dimensions, seed 2026110154 |
| Short screen | 5,000 discarded, 15,000 observed; 46 positive nonperiodic flags | Same protocol; 41 flags |
| All flagged follow-ups | 30,000 discarded, 500,000 observed; all 46 leading-positive in every window, no period through 512 detected | Same protocol; all 41 leading-positive, no period through 512 detected |
| Inspected preview pages | Four, covering every flag | Four, covering every flag |
| Longer selected cases | 109, 185, 293, 403, 414, 487, 570, 797 | 6, 7, 102, 118, 212, 673 |
| Long observations and dense portraits | Each: 100,000 discarded, five million observed; separate two-million-state portrait | Same protocol |
| Precision and 33-start controls | Cases 109, 185, 414 | Cases 6, 7 |
| Reference comparison | All eight selected cases against 108 frozen clouds | All six selected cases against the same clouds |
| Further horizon | None | Cases 6 and 7: fresh 100,000 transient and 20 million observed updates |

All fourteen selected long runs retain positive leading windows. All ten
70-/100-digit precision controls retain positive leading windows over
50,000 observed updates after 10,000 discarded steps. These controls use
exact binary64 input lifts, squared-distance Gaussian weights, nested
reflections and finite separation, without the production Jacobian.

The higher-order 33-start cohorts have 10, 14 and 11 positive nonperiodic
outcomes respectively. They also contain 22, 19 and 22 finite periodic
captures; Case 109 has one additional nonperiodic outcome with a nonpositive
leading window. Both nonuniform candidates have 33 positive nonperiodic
outcomes. These finite cohorts do not determine basin volume or uniqueness.

The nonuniform 20-million-step mean spectra are (0.308678, -0.139961,
-0.840244) and (0.357078, -0.041303, -0.405067). Every leading block is
positive. The minimum recorded normalization radii are 0.0386751 and
0.00342039. Relative to the original rendering coordinate box enlarged
by 0.02, Case 6 has zero recorded exceedances and Case 7 has two. The
threshold was chosen after the horizon calculation, before the census,
using the earlier inspection protocol. This is finite-window accounting,
not a trapping proof or invariant probability.

## Geometry and prior art

All eight preview pages and fourteen dense portraits were inspected.
Images use actual discrete states, orthogonal cameras and equal spatial
units. There is no smoothing, interpolation or symmetry replication.
All were drawn on **2026-11-01**, as recorded in the drawing-date registry.
The [current paired portrait](../../docs/assets/portraits/product-irregular.png)
and [24-second rotation](../../docs/assets/portraits/product-irregular.webp)
show the selected nonuniform cases; the article identifies their finite
observation windows. The regular higher-order images remain exploratory.
The two visual-review reports retain every image path, hash and disposition.

Small negative covariance fractions in some almost planar or linear short
records are floating eigensolver roundoff, not negative mathematical
variances. Neither the thin curves nor the folded shells are certified
surfaces, tori, fractals or closed cycles.

The comparison uses similarity normalization and affine whitening, 108
fixed reference regimes, three fingerprint shortlists, 48 starts for local
alignment and repeated subsample rankings. It does not exhaust known models
or their parameters, globally optimize shape distance, or test nonlinear
conjugacy. The [literature review](literature-review.json) records inspected
primary equations and toroidal figures. Posterior means, softmax responses,
reflection feedback, toroidal motifs and chaos in iterative algorithms all
have precedents. No historical uniqueness or strange-attractor theorem is
claimed. The earlier entropy and orientation proofs apply to Case 988 only.

## Inventory and verification

The packet inventories 25 completed run directories, preserving all 430
research files. Of these, 154 large arrays remain in the optional external
cache, with hashes, shapes and retained inputs in [omissions.json](omissions.json).
No large payload was replaced by a thinned array under its original name.
The explicitly named uniform horizon samples have stride 50.

The two complete numerical verifiers replay all 1,621 screen records,
87 flagged continuations, fourteen five-million-step calculations and all
fourteen two-million-state portraits. Histograms match exactly; independent
full-site formulas, finite-difference derivatives, local LAPACK QR and
log-determinant accounting are checked. Separate audits reconstruct both
Sobol designs exactly, replay all five 33-start cohorts and audit precision
tail returns. Both 20-million-step state histories are replayed by blocks,
with exact endpoint and normalization-minimum agreement. These replays share
the production recurrence; independent checks are identified by their scope.

[selected-summary.json](selected-summary.json), [validation.json](validation.json),
[accounting.json](accounting.json), [inventory.json](inventory.json) and
[manifest.json](manifest.json) give machine-readable results and provenance.
Run names distinguish the two surveys and each stage; the horizon and its
two excursion audits belong to the nonuniform survey.

## Reproduction

Use Python 3.12 with NumPy, SciPy, Numba, Matplotlib and mpmath 1.3.0.
Keep outputs outside the checkout. Reconstruct the earlier Cases 988/489
using the [original instructions](../product-alphabet/README.md), placing
their `case-<id>-states.npz` files in `D:/runs/original-dense` below.
Every output directory in these examples must be new.

```bash
python research/search_highorder_product.py --mode screen --source D:/runs/original-dense --out D:/runs/highorder-scout
python research/search_highorder_product.py --mode follow --source D:/runs/highorder-scout/screen.npz --out D:/runs/highorder-follow
python research/search_highorder_product.py --mode follow --source D:/runs/highorder-scout/screen.npz --out D:/runs/highorder-long --warm 100000 --count 5000000 --ids 109 185 293 403 414 487 570 797
python research/search_irregular_product.py --source D:/runs/original-dense --out D:/runs/irregular-scout
python research/search_highorder_product.py --mode follow --source D:/runs/irregular-scout/screen.npz --out D:/runs/irregular-follow
python research/search_highorder_product.py --mode follow --source D:/runs/irregular-scout/screen.npz --out D:/runs/irregular-long --warm 100000 --count 5000000 --ids 6 7 102 118 212 673
```

For each survey, use `research/render_product_alphabet.py` with its long
directory, all selected IDs and `--count 2000000`. Use a new external raw
`--out` directory and a new repository experimental `--images` directory.
Render complete previews with `research/render_projection_gallery.py`.
Run `research/verify_product_alphabet.py` with `--scout`, `--follow`, `--long`,
`--dense` and a new external `--out`.

For controls, run `research/cohort_product_alphabet.py` on the selected
case files, and `research/check_product_alphabet_precision.py` on the scout
`screen.npz` at both `--digits 70` and `--digits 100`. Set `--mpmath-path`
if mpmath is installed in a separate environment. Run
`research/compare_constellation_catalog.py` on all selected long case files;
the reference cache is documented in the [storage guide](../../docs/data-storage.md).

Use `research/audit_highorder_product.py --design highorder` or
`--design irregular`, supplying `--scout`, `--controls` (original dense
directory), `--cohort`, `--catalog`, `--long`, both `--precision` directories
and a new `--out`. Exact executed options and source snapshots are retained
in every run.

For the nonuniform horizon, run `research/extend_product_horizon.py` on its
dense directory with `--ids 6 7 --count 20000000 --block 2000000`; supply a
selection note appropriate to the replay. For each case, run
`research/inspect_product_excursions.py` with that horizon, its dense
`--reference`, `--case`, `--pad 0.02` and a new external `--out`.

# Product posterior maps: long horizon and certified finite excursions

This frozen extension follows the already published Cases 988 and 489 for
100 million updates each. It adds three exact finite escapes from the
original Case 988 portrait envelope, a certified normalization singularity,
and a full-range diagnostic. The [article](../../docs/product-alphabet/manuscript.md)
states the proofs and their limitations; earlier
[entropy](../product-symbolic/README.md) and
[orientation](../product-orientation/README.md) packets remain unchanged.

## Results and scope

| Evidence | Result |
| --- | --- |
| Case 988, fifty two-million-step QR blocks | Mean spectrum (0.251237, −0.280582, −0.440516); every leading block positive |
| Case 489, same protocol | Mean spectrum (0.337895, −0.406882, −0.787966); every leading block positive |
| Case 988 enlarged portrait envelope | 56 recorded states outside; threshold chosen after inspecting the extended bounds |
| Case 489 analogous envelope | No recorded states outside over this finite window |
| Three exact short trajectories, with open initial neighborhoods | Leave the Case 988 envelope at update 33; two interval formulas and exact rational audit agree |
| Zero of the normalization vector | Unique in the supplied tiny box, with nonsingular derivative; the map itself is undefined there |
| Three completed finer-grid proposals | All fail sampled closure; no new trapping certificate |

The threshold is the old two-million-state coordinate box enlarged by 0.02.
Its exact binary64 endpoints are stored in the two census reports. The
Case 988 frequency 56/100,000,000 is a finite-record statistic, not an
invariant probability or crisis-scaling measurement. The finite escape
certificates do not shadow the entire long orbit or identify its limit set.
The singularity is not proved to belong to that limit set. Historical
uniqueness, asymptotic chaotic attraction and equality with the certified
positive-entropy subsystem remain open.

## Inventory

All fifteen run directories are inventoried, including two unsuccessful
processes and the earlier diagnostic version. Large arrays are optional
external cache payloads; [omissions.json](omissions.json) records their hashes,
shapes and retained compact inputs. Runtime bytecode caches are separately
listed in [inventory.json](inventory.json), left outside the packet, and are
not research evidence. No numerical record was replaced by a thinned record
under its original name.

| Runs | Contents |
| --- | --- |
| `product-horizon01` | Exact inputs, fifty incoming/outgoing state and frame records, all block spectra, final-tail diagnostics, source snapshots |
| `product-excursions988-audit01`, `product-excursions489-audit01` | Every block endpoint/minimum replayed; complete threshold exceedances; explicitly named uniform stride-50 samples; fifty short minimum-centered segments per case |
| `product-excursions988-primary01`, `product-excursions988-independent01` | Forty-step interval enclosures from three selected restart boxes; factorized 140-digit and full-eight-site 180-digit formulas |
| `product-excursions988-rational01` | Exact rational endpoint audit, cross-formula intersections, six exact replay comparisons |
| `product-singularity988-primary01`, `product-singularity988-independent02` | 160-/200-digit Krawczyk certificates for a zero of 2m−x and nonsingular derivative |
| `product-singularity988-independent01` | Passed certificate written, then a final console conversion raised ValueError; corrected accepted run is independent02 |
| `product-nonuniform988-proposal01` | Initial Numba compilation failure; no successful proposal claimed |
| `product-nonuniform988-proposal02`, `product-nonuniform988-proposal03` | Monotone coordinate grids with spacing 0.02/0.01 and halfwidth 0.7/1; both enclosures leave the proposed domain |
| `product-trap988-proposal04` | Uniform spacing 1/2048; sampled enclosure leaves the domain after 38,229,697 occupied cells |
| `product-excursions988-render01`, `product-excursions988-render02` | Initial and corrected diagnostic; version 02 places cumulative minima at completed block ends and removes an overlapping label |

The three selected minimum blocks are 33, 27 and 25, in increasing order of
their minimum normalization radius. Their initial boxes have radii 1e−80
and 2e−80 in the two implementations. All three escape at update 33; the
certified margins exceed 0.46319, 0.40401 and 0.18140. The root box for the
normalization singularity has radius 1e−50 in the primary calculation and
2e−50 in the full-site check. The latter Krawczyk image lies inside the
primary root box. Both formulations share mpmath's interval backend and a
Newton center proposal; independence is at the expression and derivative
level, not the interval-library level.

[selected-summary.json](selected-summary.json),
[validation.json](validation.json), [accounting.json](accounting.json) and
[manifest.json](manifest.json) give compact results and byte-level provenance.
The [literature receipt](literature-review.json) records the additional
attention/normalization comparison, inspected figures and abstract-only
crisis precedent. Third-party screenshots are not redistributed.

## Reproduction

Use Python 3.12, NumPy, Numba, SciPy, Matplotlib and mpmath 1.3.0. The exact
executed environment versions appear in the validation record. Keep generated
arrays and replay outputs outside the checkout. Reconstruct the original
`case-988-states.npz` and `case-489-states.npz` using the
[original campaign instructions](../product-alphabet/README.md); expected
hashes are in `product-horizon01/protocol.json`. The following examples assume
those files are in `D:/runs/product-dense` and use new output directories.

```bash
python research/extend_product_horizon.py --source D:/runs/product-dense --out D:/runs/product-horizon --count 100000000 --block 2000000
python research/inspect_product_excursions.py --horizon D:/runs/product-horizon --reference D:/runs/product-dense/case-988-states.npz --case 988 --out D:/runs/excursions988
python research/inspect_product_excursions.py --horizon D:/runs/product-horizon --reference D:/runs/product-dense/case-489-states.npz --case 489 --out D:/runs/excursions489
python research/certify_product_excursions.py --source D:/runs/excursions988 --out D:/runs/escapes-primary --engine factorized --digits 140 --radius 1e-80
python research/certify_product_excursions.py --source D:/runs/excursions988 --out D:/runs/escapes-full --engine full --digits 180 --radius 2e-80
python research/certify_product_singularity.py --input results/product-symbolic/runs/product-coverings988-certificate04/input.json --out D:/runs/singularity-primary --engine factorized --digits 160 --radius 1e-50
python research/certify_product_singularity.py --input results/product-symbolic/runs/product-coverings988-certificate04/input.json --out D:/runs/singularity-full --engine full --digits 200 --radius 2e-50
```

For the exact endpoint audit and six payload replays, set `--mpmath-path` to
the environment's site-packages directory. The audit does not rerun the
hundred-million-step calculations.

```bash
python research/audit_product_excursions.py --source D:/runs/excursions988 --out D:/runs/escape-audit --replay D:/runs/escape-replays --mpmath-path D:/runs/env/Lib/site-packages --singularities D:/runs/singularity-primary D:/runs/singularity-full --excursions D:/runs/escapes-primary D:/runs/escapes-full
```

For the diagnostic, use `research/render_product_excursions.py` with
`--source D:/runs/excursions988`, `--camera` pointing to the original
`case-988.json`, a new external `--out`, and an image path inside the repository.
Its retained `figure-data.npz` contains the complete histogram, all projected
exceedances, event times, block minima and camera. Histograms and events were
recomputed exactly before freezing. The blue density uses every fiftieth
state; all 56 threshold exceedances are overlaid separately. No smoothing,
nonlinear coordinate transformation, clipped plotted points or invented
connecting trajectory is used.

Grid proposal parameters and full progress histories are retained in each
run's report. Their inflated sampled corner residuals are not rigorous
Taylor remainder bounds. Reproducing their failure does not prove the
nonexistence of every possible trapping region.

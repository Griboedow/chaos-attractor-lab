# Circle / finite-point RAAR campaign

Numerical exploration on 2026-11-01. See the [scientific report](../../docs/raar-constellations.md) for the model, selection, controls and prior-art limits.

`inventory.json` accounts for all 172 files in eight completed calculation directories. It includes intermediate scout blocks and the partial comparison report. `omissions.json` records 104 large optional arrays, their exact hashes and retained compact keys. Original full arrays remain in their external run directories and the hash-addressed cache; none is replaced by a thinned array under the same name. The compact archives preserve exact inputs and small metrics, not full trajectories. See [storage](../../docs/data-storage.md).

| Run | Purpose |
|---|---|
| `raar-constellation-scout01` | All 2,048 geometric cases; 72 flags; independent algebra and finite-difference checks |
| `raar-constellation-follow01` | Every flag, fresh 500,000-step follow-up; periodic captures retained |
| `raar-constellation-long01` | Six visually selected cases, five million measured iterations each |
| `raar-constellation-precision01` | Cases 1504/1725, scalar nested reflectors at 70/100 digits |
| `raar-constellation-cohort01` | Six cases, 33 initial states each |
| `raar-constellation-render01` | Six full two-million-state clouds, histograms and fixed cameras |
| `raar-constellation-catalog01` | Similarity and affine-normalized comparisons with 108 reference regimes |
| `raar-constellation-branches01` | Three-site reduction and 300-value beta sweeps |

The first 72-case portrait attempt failed before producing images because its camera matrix was not transposed. The failure record and executed source are retained in `presentation-failure/`; corrected previews are in `docs/assets/experimental/constellation-gallery02`. All six dense portraits were inspected, with Cases 1504/1725 selected for the main-page rotation. The aesthetic selection is declared and does not alter earlier numerical outcomes.

## Recompute

Use Python 3.12, NumPy, SciPy, Numba, Matplotlib, Pillow and mpmath; exact observed versions are in `environment.json`. Set Numba and Matplotlib caches outside a small system disk. Run from the repository root. Every output directory must be new. The commands below use a separate D: directory; change that prefix on other platforms.

```text
python research/search_raar_constellations.py --out D:/Codex/new-constellation-scout --power 9
python research/follow_raar_constellations.py --source D:/Codex/new-constellation-scout/screen.npz --out D:/Codex/new-constellation-follow
python research/extend_raar_constellations.py --source D:/Codex/new-constellation-scout/screen.npz --out D:/Codex/new-constellation-long --ids 124 611 896 1504 1725 1873 --count 5000000
python research/check_raar_constellation_precision.py --source D:/Codex/new-constellation-scout/screen.npz --out D:/Codex/new-constellation-precision --ids 1504 1725 --digits 70 100
python research/cohort_raar_constellations.py --source D:/Codex/new-constellation-long --out D:/Codex/new-constellation-cohort --ids 124 611 896 1504 1725 1873
python research/render_raar_constellations.py --source D:/Codex/new-constellation-long --out D:/Codex/new-constellation-render --images docs/assets/experimental/new-constellation-render --ids 124 611 896 1504 1725 1873
python research/control_constellation_branches.py --source D:/Codex/new-constellation-render --out D:/Codex/new-constellation-branches --ids 1504 1725
python research/compare_constellation_catalog.py --source D:/Codex/new-constellation-long --out D:/Codex/new-constellation-catalog --ids 124 611 896 1504 1725 1873
```

The catalog comparison also needs the frozen reference arrays, recoverable with `python tools/data_cache.py fetch --group attractor-novelty`. Recomputed chaotic trajectories can differ across library/compiler environments; original cache hashes identify the historical run. The current scripts and executed snapshots are both retained. None of the optional arrays in this new campaign has a historical Git download.

With the original cache present, `python research/verify_constellation_campaign.py --out D:/Codex/new-constellation-verification` checks the inventory, repeats all six five-million-step runs and all six dense trajectories, and reconstructs every dense histogram. This is a numerical replay audit, not a computer-assisted proof of chaos. The observed verification result is in `verification/`.

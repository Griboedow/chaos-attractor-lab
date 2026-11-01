# From a tiny segment to a serrated arch

[Read the six-page note](manuscript.pdf) · [Watch the evolution](../../docs/assets/raar-cycle-unfolding.gif) · [Illustrated report](../../docs/raar-folding.md) · [Base manuscript](../raar-circles/manuscript.pdf) · [Cycle companion](../raar-cycle-structure/manuscript.pdf)

The note contains four figures, six complete finite-evolution controls, a degree-eight inverse reduction, complete preimage counts for 25 exact targets, an independent interval verification and a carefully limited conjugacy exclusion. The appendix checks the seven time-residue samples; the updated prior-art discussion credits an explicit attractor portrait in a 2001 preprint. It is a research note, not a peer-reviewed publication or a declaration of historical priority.

The supported statements are separated from unresolved ones: tangent segments are not certified unstable manifolds; finite sample agreement is not asymptotic attraction; full-map noninjectivity is not noninjectivity on the attractor; no exact literature match is not a proof of novelty.

## Reproduce

Run from the repository root, using a fresh output directory for each command:

```sh
python research/verify_raar_folding.py
python research/verify_raar_phase_review.py
python research/audit_raar_folding_manuscript.py
python research/trace_raar_cycle_seeds.py --out outputs/raar-seeds-replay
python research/raar_seed_hires_histograms.py --source outputs/raar-seeds-replay --out outputs/raar-seeds-hires-replay
python research/certify_raar_inverse.py --input results/raar-folding/inverse-proof/input.json --out outputs/raar-inverse-replay
python research/check_raar_preimages_interval.py --input results/raar-folding/inverse-independent/input.json --sturm outputs/raar-inverse-replay/report.json --out outputs/raar-inverse-independent-replay
python research/build_raar_folding_manuscript.py --output outputs/raar-folding-paper-replay/manuscript.pdf
```

Python 3.12.14; numerical dependencies: NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1, Matplotlib 3.11.2 and Pillow 13.0.0. Interval certificates use mpmath 1.3.0 at 100 and 200 decimal digits. Paper generation uses ReportLab 4.5.1 and pypdf 6.21.4. Interval proofs trust the arithmetic library; they are not proof-assistant formalizations.

The evidence [manifest](../../results/raar-folding/manifest.json) covers all six experiments, both numerical inverse searches, all certificate revisions, all rendering revisions and provenance. The final rendering supersedes earlier resampling wording. This folder contains editable prose, exact figure copies, a build record, a page-review record and its own manifest. Every page was inspected after Poppler rendering at 110 dpi; scientific figure raster layers are embedded at 300 dpi.

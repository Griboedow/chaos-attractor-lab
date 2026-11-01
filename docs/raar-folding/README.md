# From a tiny segment to a serrated arch

[Read the note](manuscript.md) · [Watch the evolution](../assets/raar-cycle-unfolding.gif) · [Illustrated report](../raar-folding.md) · [Base manuscript](../raar-circles/manuscript.md) · [Cycle companion](../raar-cycle-structure/manuscript.md)

The note contains five figures, six complete finite-evolution controls, a degree-eight inverse reduction, complete preimage counts for 25 exact targets, an independent interval verification and a carefully limited conjugacy exclusion. Two appendices check the seven time-residue samples and certify a smooth fold near the finite reference cloud. The prior-art discussion credits an explicit attractor portrait in a 2001 preprint. It is a research note, not a peer-reviewed publication or a declaration of historical priority.

The supported statements are separated from unresolved ones: tangent segments are not certified unstable manifolds; finite sample agreement is not asymptotic attraction; full-map noninjectivity is not noninjectivity on the attractor; no exact literature match is not a proof of novelty.

## Reproduce

Run from the repository root, using a fresh output directory for each command:

```sh
python research/verify_raar_folding.py
python research/verify_raar_phase_review.py
python research/verify_raar_fold_certificate.py
python research/trace_raar_cycle_seeds.py --out outputs/raar-seeds-replay
python research/raar_seed_hires_histograms.py --source outputs/raar-seeds-replay --out outputs/raar-seeds-hires-replay
python research/certify_raar_inverse.py --input results/raar-folding/inverse-proof/input.json --out outputs/raar-inverse-replay
python research/check_raar_preimages_interval.py --input results/raar-folding/inverse-independent/input.json --sturm outputs/raar-inverse-replay/report.json --out outputs/raar-inverse-independent-replay
```

Python 3.12.14; numerical dependencies: NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1, Matplotlib 3.11.2 and Pillow 13.0.0. Interval certificates use mpmath 1.3.0 at 100 and 200 decimal digits. Interval proofs trust the arithmetic library; they are not proof-assistant formalizations.

The evidence [manifest](../../results/raar-folding/manifest.json) covers all six experiments, both numerical inverse searches, all certificate revisions, all rendering revisions and provenance. Separate packets retain the [phase diagnostic](../../results/raar-phase-classes), [extended review](../../results/raar-priority-extension) and [smooth-fold certificate](../../results/raar-fold-certificate).

The complete text is edited directly in [manuscript.md](manuscript.md). Run `python tools/check_manuscripts.py` after editing. The original build, review and manifest records describe the historical publication edition; see [editing and archive instructions](../manuscript-editing.md). Frozen numerical verifiers may require the explicitly restored files described there.

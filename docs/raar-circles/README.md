# Serrated arches and layered loops in relaxed reflections between spatial circles

[Read the illustrated manuscript](manuscript.md) · [Rotate the actual states](../assets/raar-circles-rotation.gif) · [Public report](../raar-circles.md) · [Exact structure](../raar-structure.md) · [Priority audit](../raar-priority.md)

Six figures, editable prose and reproducible calculations. This is a research manuscript, not a peer-reviewed publication or a priority certificate.

## Supported claims

- Two selected spatial-circle configurations of the existing RAAR algorithm produce intricate finite state clouds, reproduced with independent high-precision arithmetic and the original ProxToolbox evaluator.
- Five-million-step spectra support numerical strange-attractor candidates. R820 has two positive mean rates but one negative short second-rate window in the original-library check.
- An elementary absorbing bound applies to every defined orbit. Interval Sturm calculations completely classify fixed points at six exact settings: five have none; R1679 has exactly two, both unstable.
- All six input-circle pairs are disjoint. The R1679 pair has certified linking number −1; this is not a topology claim about its attracting cloud.
- The packet retains every declared screen/follow-up result, periodic neighbors, all 53 gallery clouds, 384 cohort experiments, 131,072 initial-volume trajectories and 771 parameter settings.

RAAR, planar-circle examples, projection-search chaos and general fixed-point relations are prior work. No new algorithm, new general chaos mechanism, asymptotic strange attractor or historically unique attractor class is claimed.

## Reproduction

From the repository root, using Python 3.12:

```text
python research/verify_raar_circles.py --root results/raar-circles --assets docs/assets --output NEW_PACKET_REVIEW.json
```

The numerical environment uses NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1, Matplotlib 3.11.2 and Pillow 13.0.0. The certificate environment uses mpmath 1.3.0 and SymPy 1.14.0.

```text
python research/certify_raar_fixed_points.py --input results/raar-circles/fixed-proof/input.json --out NEW_CERTIFICATE_DIRECTORY
python research/audit_raar_algebra.py --input results/raar-circles/algebra/input.json --out NEW_ALGEBRA_DIRECTORY
python research/search_raar_circles.py --out NEW_SCOUT_DIRECTORY
python research/follow_raar_circles.py --source results/raar-circles/scout/screen.npz --out NEW_FOLLOW_DIRECTORY
python research/follow_raar_circles.py --source results/raar-circles/scout/screen.npz --out NEW_LONG_DIRECTORY --ids 75 820 1051 1311 1483 1679 --steps 5000000
python research/check_raar_precision.py --cases results/raar-circles/precision/input.json --out NEW_PRECISION_DIRECTORY
python research/check_raar_upstream.py --source results/raar-circles/scout/screen.npz --upstream results/raar-circles/upstream-source --out NEW_LIBRARY_DIRECTORY
```

Use new output directories: calculation scripts refuse to overwrite existing runs. Source snapshots and complete protocols accompany the arrays. Exact parameter values are stored in the arrays and certificate input; printed decimals are shortened.

The [packet manifest](../../results/raar-circles/manifest.json) covers the evidence. Original third-party source is retained with its license; publication PDFs are linked rather than redistributed. The Borwein et al. article has a 2017 online DOI record and a 2018 journal-volume citation, 26:385–403; the manuscript uses the latter.

The complete text is edited directly in [manuscript.md](manuscript.md). Run `python tools/check_manuscripts.py` after editing. The original build, review and manifest records describe the historical publication edition; see [editing and archive instructions](../manuscript-editing.md). Frozen numerical verifiers may require the explicitly restored files described there.

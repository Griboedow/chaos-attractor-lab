# Periodic structure behind a serrated spatial RAAR cloud

[Read the six-page companion](manuscript.pdf) · [Illustrated report](../../docs/raar-cycle-structure.md) · [Base manuscript](../raar-circles/manuscript.pdf) · [Rotate the actual cloud](../../docs/assets/raar-circles-rotation.gif) · [Updated priority audit](../../docs/raar-priority.md)

Three figures, editable prose, independent interval calculations and a complete fine parameter slice. This is a research companion note, not a peer-reviewed publication or a priority certificate.

The supported results are four attracting cycles of least periods 7/14/28/56 at different settings, one unstable least-seven cycle at the base arch, and four brackets containing simple −1 multiplier events. At each event the other two multipliers lie inside the unit disk; its local branch is attracting at the left endpoint and unstable at the right. Neither nondegenerate flip conditions, an infinite cascade, asymptotic chaos nor historical uniqueness is proved.

The base seven-cycle has one unstable real multiplier and two stable real multipliers, with a negative full-period determinant. The standard orientation invariant excludes a smooth conjugacy near that cycle to a fixed-time ODE solution map on R³. This does not establish inequivalence to every discrete attractor or classify all return-map and higher-dimensional constructions.

## Reproduce and audit

Run from the repository root; use new output directories.

```text
python research/verify_raar_cycles.py --output NEW_PACKET_REVIEW.json
python research/build_raar_cycle_manuscript.py --output NEW_BUILD_DIRECTORY/manuscript.pdf
python research/audit_raar_cycle_manuscript.py

python research/certify_raar_cycles.py --input results/raar-cycle-structure/cycle-proof/input.json --out NEW_CYCLE_CERTIFICATES
python research/certify_raar_multiplier_events.py --input results/raar-cycle-structure/multiplier-proof/input.json --out NEW_EVENT_CERTIFICATES
python research/audit_raar_cycle_structure.py --source results/raar-cycle-structure/cycle-proof/base-unstable-period-7.json --out NEW_STRUCTURE_AUDIT

python research/continue_raar_cycles.py --source results/raar-circles/relaxation/scan.npz --out NEW_CONTINUATION
python research/scan_raar_cascade.py --source results/raar-circles/scout/screen.npz --out NEW_FINE_SLICE
python research/build_raar_cycle_figures.py --root results/raar-cycle-structure --out NEW_FIGURES
```

The interval environment uses Python 3.12 and mpmath 1.3.0, with 200 digits for isolated cycles and 600 digits for uniform event brackets. Interval results trust this arithmetic library; they are not proof-assistant formalizations. Numerical/rendering dependencies are NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1, Matplotlib 3.11.2, ReportLab 4.5.1 and pypdf 6.21.4. Every PDF page is reviewed after Poppler rendering at 110 dpi; rasterized scientific figure layers are embedded at 300 dpi.

The [packet manifest](../../results/raar-cycle-structure/manifest.json) identifies all evidence, failed/earlier attempts and source snapshots. The original spatial-circle packet is unchanged. This folder's manifest covers the PDF, editable content, three figures, build record and visual review. The bounded prior-art review explicitly credits existing RAAR, planar circles, spatial DR nonconvergence, projection-search chaos, DR/ADMM chaos conjectures, period doubling and interval verification.

# Gaussian-posterior attracting dynamics: symbolic and periodic certificates

The [article](../../docs/soft-symbolic-dynamics/manuscript.md) proves positive
entropy inside the [previously certified attracting set](../../docs/soft-constellation-trapping.md),
locates periodic points with different unstable dimensions, and derives two
restricted C1 equivalence obstructions. Historical uniqueness is unresolved.

## Retained evidence

- [Primary covering certificate](runs/soft-covering-certificate01/report.json):
  seven one-step crossings, 21 complete volume/face checks, 4,737 leaves.
- [Independent automatic-derivative certificate](runs/soft-covering-ad01/report.json):
  direct Gaussian-distance map, 80-digit intervals, 9,918 leaves.
- [Primary exact audit](runs/soft-covering-audit01/report.json) and
  [second exact audit](runs/soft-covering-ad-audit01/report.json): complete
  dyadic tilings, strict inequalities, exact matrix inverses, disjoint
  windows, and membership of N0 in the existing trapping region.
- [Periodic-point certificates](runs/soft-cycle-certificates01/report.json)
  and [independent AD certificates](runs/soft-cycle-ad01/report.json):
  200-digit unique-root inclusions, least periods, all three multiplier
  brackets and full phase membership in the region.
- [Cross-packet verification](runs/soft-symbolic-verification01/report.json):
  identities, source hashes, strict saved inequalities, distinct cycles and
  final display coordinates.
- [Bounded prior-art review](runs/soft-symbolic-priority01/report.json) and
  [visual review](runs/soft-symbolic-priority01/presentation-review.json).

The packet inventories all 199 original files in 13 runs. The single large
return-scout array, 13,282,450 bytes, remains outside Git with its hash and
compact inputs retained. See `inventory.json`, `omissions.json`, and the
[storage policy](../../docs/data-storage.md). All 640 numerical root attempts
and all 39 window proposals are retained, including failed sampled proposals.
The first rendering revision has an overlapping explanatory note; its source
and image remain archived, and the corrected second revision is published.

The complete window input and all cycle seed files are small and retained.
Checking the mathematical certificates requires no long trajectory arrays.
The earlier trapping and numerical packets remain unchanged.

## Reproduce the principal certificates

Use Python 3.12, NumPy 2.5.4, Numba 0.68.1 and mpmath 1.3.0, matching the
executed environment. SciPy 1.18.1 is needed for discovery, and Matplotlib
3.11.2 for rendering. Every output directory below must be new. An optional
`--mpmath-path` appends a separate library directory without overriding the
NumPy environment.

```text
python research/certify_soft_coverings.py --input results/soft-symbolic-dynamics/runs/soft-covering-proposals01/windows-1.npz --trap results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --out D:/Codex/new-soft-covering
python research/audit_soft_coverings.py --source D:/Codex/new-soft-covering --trap results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --out D:/Codex/new-soft-covering-audit
python research/check_soft_symbolic_ad.py --source D:/Codex/new-soft-covering --out D:/Codex/new-soft-covering-ad
python research/audit_soft_coverings.py --source D:/Codex/new-soft-covering-ad --trap results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --out D:/Codex/new-soft-covering-ad-audit
python research/certify_soft_cycles.py --source results/soft-symbolic-dynamics/runs/soft-cycles01 --trap results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --out D:/Codex/new-soft-cycles
python research/check_soft_cycles_ad.py --source D:/Codex/new-soft-cycles --trap results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --out D:/Codex/new-soft-cycles-ad
```

Saved partitions can be audited directly by replacing the new source paths
with the corresponding retained run directories. The second covering engine
recomputes every first-engine leaf and may subdivide further. Its `examined`
field counts nodes in the complete final binary tree for the common audit;
`recomputation_calls` records actual interval calls in the second engine.
These different counts are not run-time comparisons.

The root certificates use radius `1e-100` in the primary implementation and
`1e-90` in the independently unique-root control. Three disjoint real
multiplier brackets exhaust the cubic's roots. Exact membership and
periodicity put all certified phases in the maximal invariant set K.

## Reproduce discovery and figures

Recover or recompute Case 156 from the
[numerical campaign](../soft-constellations/README.md) for discovery and the
two-million-state rendering. Discovery is not a proof prerequisite once
the fixed windows and cycle seeds are supplied.

```text
python research/find_soft_cycles.py --source D:/Codex/soft-long/case-156.npz --trap results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --out D:/Codex/new-soft-cycle-scout
python research/search_soft_returns.py --source D:/Codex/new-soft-cycle-scout/cycle-0.npz --out D:/Codex/new-soft-return-scout
python research/propose_soft_coverings.py --source D:/Codex/new-soft-return-scout --out D:/Codex/new-soft-window-proposals
python research/render_soft_symbolic.py --states D:/Codex/soft-render/case-156-states.npz --view D:/Codex/soft-render/case-156.json --cycles D:/Codex/new-soft-cycles --out D:/Codex/new-soft-symbolic-render --figures D:/Codex/new-soft-symbolic-figures
python research/verify_soft_symbolic.py --base results/soft-symbolic-dynamics/runs --trap results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --out D:/Codex/new-soft-symbolic-cross-check
```

The cross-check uses the archived run names and published figure bytes.
The source snapshots in each run identify the actual executed versions.
The tangent scout is not a validated unstable manifold; small return
residuals and sampled window tests are not the covering proof.

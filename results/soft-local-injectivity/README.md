# Connectedness and local noninjectivity: Case 156

This packet checks the finite statements supporting the topological result
in [Section 6.1 of the article](../../docs/soft-symbolic-dynamics/manuscript.md#61-an-obstruction-under-continuous-changes-of-ambient-coordinates).
It depends on the frozen [trapping packet](../soft-trapping/README.md) and
[periodic-point certificates](../soft-symbolic-dynamics/README.md).

All 35,784 trapping cubes form one component under full-face adjacency.
The complete spanning tree has 35,783 edges; a separate union-find algorithm
finds the same component and 88,171 undirected face adjacencies. The original
cube list and grid are retained in the trapping packet. This packet retains
the entire tree, not a sample.

Exact rational interval determinants of both independently computed
periodic return matrices prove det DF(P) > 0 and det DF⁵(A) < 0. A separate
audit confirms the signs using the three certified multiplier brackets.
The fixed point and every phase of the five-cycle belong to the attracting
set K. The negative return determinant forces a negative one-step
determinant at at least one phase. Approximate values are display aids;
the reports contain the exact rational enclosures.

Since Q is connected, its nested continuous images have connected
intersection K. Opposite orientation signs on K contradict local
injectivity everywhere on K, by invariance of domain and the locally
constant orientation of a local homeomorphism. This final deduction is a
written topological proof; it is not machine formalized.

The result excludes ambient neighborhood conjugacy to a locally injective
map. It does **not** certify two distinct points of K with the same image,
a fold normal form, transitivity, or historical uniqueness.

## Retained evidence

- [Primary certificate](runs/soft-local-injectivity01/report.json), executed
  source and complete compressed spanning tree.
- [Independent audit](runs/soft-local-injectivity-verification01/report.json)
  and its executed source.
- [Inventory](inventory.json), [manifest](manifest.json) and
  [bounded literature review](literature-review.json).

There are no omitted arrays in these two runs. Input hashes bind the exact
trapping grid and both original periodic-point reports. The earlier packets
are dependencies and remain unchanged.

## Reproduction

Use Python 3.12 and NumPy 2.5.4. The connectivity and orientation checks use
integer addresses and Python's exact rational arithmetic; they do not
depend on a long trajectory or a new numerical root search. Choose fresh
output directories:

```text
python research/certify_soft_local_injectivity.py --trap results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --cycles results/soft-symbolic-dynamics/runs/soft-cycle-certificates01/report.json --control results/soft-symbolic-dynamics/runs/soft-cycle-ad01/report.json --out D:/Codex/new-soft-local-injectivity
python research/verify_soft_local_injectivity.py --source D:/Codex/new-soft-local-injectivity --trap results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --cycles results/soft-symbolic-dynamics/runs/soft-cycle-certificates01/report.json --control results/soft-symbolic-dynamics/runs/soft-cycle-ad01/report.json --out D:/Codex/new-soft-local-injectivity-audit
```

The full trapping and root-existence calculations are reproduced through
their respective packets. These shorter checks do not replace them.

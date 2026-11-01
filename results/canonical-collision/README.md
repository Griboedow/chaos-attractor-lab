# Canonical cube: intrinsic collision and homoclinic return

The [article, Section 5](../../docs/canonical-sphere/manuscript.md#5-certified-attraction-merged-histories-and-a-homoclinic-return)
gives the theorem and graph-transform argument. The exact stored Case 653
has a strict cubic trapping region. Its maximal attracting set contains
distinct points with the same image and a nontrivial full orbit returning
to the same symmetric fixed point in both time directions. Noninjectivity
excludes conjugacy of this whole set with an injective map.

The theorem does not identify this maximal attracting set with the numerical
portrait, or establish asymptotic hyperchaos, transitivity, a horseshoe,
positive entropy or historical uniqueness. Symmetric noninvertible maps
and homoclinic geometry are established mechanisms; see the
[priority review](../../docs/soft-constellation-priority.md#intrinsic-collisions-and-homoclinic-geometry)
and [literature receipt](literature-review.json).

## Accepted certificates

| Component | Primary record | Independent formula check |
| --- | --- | --- |
| Strict trapping, fixed point, local unstable graph, collision | [250-digit certificate](runs/canonical-collision04/report.json) | [300-digit eight-site automatic derivatives](runs/canonical-collision-verification02/report.json) |
| Homoclinic crossing and stable return | [421-step crossing; 61-step return](runs/canonical-homoclinic04/report.json) | [Eight-site verification; 121-step return to a smaller square](runs/canonical-homoclinic-verification02/report.json) |
| Diagnostic figure | [Final rendering record](runs/canonical-collision-figure03/report.json) | Visual review of the rendered figure; it is not itself a certificate |

All input binary64 values are interpreted as exact real numbers. The cube
factorization uses the exact half-sum and half-difference of the two stored
site endpoints. Fixed-point boxes include binary endpoint tuples, and the
coordinate matrix and inverse are exact rational numbers. The unstable
graph is enclosed rather than replaced by its tangent line.

The separate collision verifier implements eight unshifted Gaussian weights,
nested reflections and automatic differentiation without importing the
production factorization or Jacobian. It uses 300 rather than 250 decimal
digits and an 8,192-cell rather than 4,096-cell scalar trapping partition.
The homoclinic check shares the propagation driver and nominal centers with
the primary calculation, but replaces every enclosing map/Jacobian
evaluation by eight-site automatic derivatives. It also uses a smaller
stable square and a longer validated return. Both use mpmath 1.3.0 interval
arithmetic, so this is formula independence, not arithmetic-library independence.

## Complete history and unsuccessful attempts

Sixteen run directories contain **51 inventoried files**, all retained.
No new large trajectory arrays are needed. The original numerical campaign
and its full source input remain frozen in [canonical-sphere](../canonical-sphere/README.md).

The 32-phase [search](runs/canonical-unstable-search01/report.json) is only
a numerical guide. Its straight initial line is not asserted to be the
unstable graph. Runs `canonical-collision02` and `03` are earlier successful
iterations; run `04` uses exact point conversion and exports binary root
endpoints. The first verifier was superseded by verification `02`, which
also checks inclusion of the entire primary root box in its independently
validated uniqueness box.

The failure receipts retain the actual scope of each unsuccessful attempt:

- `canonical-collision01`: a natural scalar quotient enclosure was too wide
  to establish strict trapping; monotonic endpoint bounds resolve it.
- `canonical-homoclinic01` and `02`: the assertion assumed the wrong
  crossing orientation. The selected crossing has positive left and
  negative right endpoint values.
- `canonical-homoclinic03`: direct Cartesian propagation of the planar
  return wrapped too widely. Centered stable coordinates resolve it.
- `canonical-homoclinic-verification01`: the eight-site derivative bound
  on the radius `1e-4` square was too wide to prove contraction. Verification
  `02` uses radius `1e-6`, enters after 121 steps, and proves contraction.

The three figure runs retain successive presentation changes. The final
asset is a diagnostic of enclosure centers; dashed segments are discrete
iteration guides. It does not draw a continuous flow or claim that the
homoclinic orbit fills the sampled density portrait.

## Reproduction

Use Python 3.12, NumPy and mpmath 1.3.0. Matplotlib is only needed for the
diagnostic. Supply the installed mpmath package directory to `--mpmath-path`.
Use fresh external output directories for every run. The exact source array
has SHA-256
`750276b784b4ba40a1f6836d119cf2b294d84e65f841d80cb33d17745e0849e5`.
It is registered as
`results/canonical-sphere/runs/canonical-sphere-scout01/screen.npz`.
Use [data storage](../../docs/data-storage.md) to resolve the cache, or follow
the original campaign's reproduction instructions. The retained compact
input has the exact parameter/site arrays but a different archive hash;
the commands below assume the full original array at `D:/Codex/input/screen.npz`.

```text
python research/search_canonical_unstable.py --source D:/Codex/input/screen.npz --out D:/Codex/new-unstable-search --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages
python research/certify_canonical_collision.py --source D:/Codex/input/screen.npz --search D:/Codex/new-unstable-search --out D:/Codex/new-collision --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages
python research/verify_canonical_collision.py --source D:/Codex/input/screen.npz --certificate D:/Codex/new-collision/report.json --out D:/Codex/new-collision-verification --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages
python research/certify_canonical_homoclinic.py --source D:/Codex/input/screen.npz --certificate D:/Codex/new-collision/report.json --out D:/Codex/new-homoclinic --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages
python research/certify_canonical_homoclinic.py --source D:/Codex/input/screen.npz --certificate D:/Codex/new-collision/report.json --out D:/Codex/new-homoclinic-verification --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages --independent
python research/render_canonical_collision.py --collision D:/Codex/new-collision/report.json --homoclinic D:/Codex/new-homoclinic/report.json --out D:/Codex/new-collision-figure --image docs/assets/diagnostics/new-canonical-collision.png
python tools/render_dates.py record docs/assets/diagnostics/new-canonical-collision.png --role diagnostic
```

Executed snapshots accompany every run. NumPy selects a convenient stable
basis, which is then treated as exact rational data and validated. A different
linear-algebra environment may select different eigenvector signs or order;
the inequalities, not bitwise identity of that choice, are the test.
`inventory.json` accounts for every file in the completed or failed runs.
`manifest.json` binds this packet's files, excluding the manifest itself.

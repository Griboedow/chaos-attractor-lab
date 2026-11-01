# Positive entropy in the four-wing product-alphabet map

The [article, Section 6](../../docs/product-alphabet/manuscript.md#6-a-certified-positive-entropy-subsystem)
proves a compact invariant subsystem C for the exact Case 988 parameters,
with **topological entropy above 0.0254 per original iteration**. Twenty-four
coverings of G = F^4 create returns of lengths one and 23. The positive root
of `rho^23 - rho^22 - 1 = 0` gives the bound `log(rho)/4`.

This is a certificate of recurrent complexity, not a trapping or attraction
certificate. Equality of C with the four-wing numerical cloud, uniform
hyperbolicity, a transverse homoclinic orbit and historical uniqueness remain
unproved. The [numerical survey](../product-alphabet/README.md) is a separate,
unchanged packet. Covering relations and irregular reflected algorithms are
known; see the [supplementary literature receipt](literature-review.json).

## Accepted records

| Component | Record |
| --- | --- |
| Exact windows, parameters and eight sites | [Rational coordinate data](runs/product-coverings988-certificate04/input.json) |
| 80-digit factorized interval calculation | [Primary certificate](runs/product-coverings988-certificate04/report.json), [exact audit](runs/product-coverings988-audit04/report.json) |
| 110-digit full-eight-site automatic derivatives | [Independent formula check](runs/product-coverings988-independent04/report.json), [exact audit](runs/product-coverings988-audit05/report.json) |
| Identity with published Case 988 inputs | [Bitwise input and dependency check](runs/product-coverings988-input-check01/report.json) |
| Diagnostic window locations and graph | [Final rendering record](runs/product-coverings988-render02/report.json) |

The primary and second calculations retain 178 and 210 accepted leaves,
respectively. Audits use exact binary interval endpoints, complete binary
subdivision trees, exact rational matrix inverses and a separating coordinate
for N0 and N1. Minimum normalized margins exceed 0.00629 and 0.01915; the
normalization radius exceeds 2.3799 throughout the validated four-step
enclosures. Both calculations use mpmath 1.3.0 intervals, so independence
concerns map formulas and differentiation, not the arithmetic library.

The theorem uses closed-chain existence from Zgliczynski and Gidea (2004).
Words formed from returns of length one and 23 obey `a[n]=a[n-1]+a[n-23]`.
Their first different choice separates the realizing states between N0 and
N1. Other windows need not be disjoint. The closure of the resulting periodic
orbits is compact and invariant in the smooth domain. Attraction is not
inferred from their proximity to sampled states.

## Complete exploratory record

Thirty runs contain **643 inventoried files**. Of these, 636 raw files are
retained, five large arrays remain in the optional cache, and two incomplete
interval-leaf logs remain external. Their hashes and locations are recorded
in [external-logs.json](external-logs.json). Neither external log is needed
to check an accepted certificate. The complete [inventory](inventory.json)
and [run accounting](accounting.json) include failures and superseded runs.

- `product-cycles988-01`: 512 numerical fits across periods 1–16, producing
  70 cycle proposals. Cycle 14, of proposed period four, guides the search.
  No interval cycle-existence or manifold assertion is imported from it.
- `product-returns988-01`, `refined01`, `prepend2`, `prepend4`: unstable-line
  searches and multiple shooting with linear endpoint constraints. These
  propose window locations, not validated homoclinic connections.
- `product-trap988-proposal01` and `02`: nonclosing grid proposals. They do
  not establish a trapping region.
- `product-coverings988-01`, `03`, `prepend2`, `prepend4`: unsuccessful
  sampled proposals. Run `02` stopped on a list/array construction error.
- `product-coverings988-scaled2` and `scaled4`: the linear program is solved
  in microunits, with tighter tolerances and explicit constraint residuals.
  All 24 sampled proposals pass; only selected windows receive interval checks.
- `product-coverings988-certificate01`: incomplete because a natural
  quotient interval was too wide. Runs `02` and `03` pass using an equivalent
  stable logistic expression. Run `04` adds exact binary endpoint exports.
- `product-coverings988-independent01`: an uncentered full-site expression
  was terminated incomplete because of interval overestimation. Run `02`
  stopped on operand dispatch between an interval and an automatic-derivative
  object. Run `03` passes after the correction; final run `04` additionally
  checks that each mean-value center belongs to its current box.
- Five audit runs and two figure runs preserve successive checks and layouts.
  The final audits add the sharper return-count entropy bound. The figure is
  a diagnostic; the original two-million-state arrays remain in the numerical
  campaign's external cache.

Large proposal reports are retained as explicitly linked evidence of the
complete selection history. Raw arrays are not replaced by smaller arrays
under their original names. [Data storage](../../docs/data-storage.md) explains
cache resolution and the separately named compact records.

## Reproduction

Use Python 3.12, NumPy, SciPy, Numba and mpmath 1.3.0. Matplotlib is needed
only for the diagnostic. Supply the directory containing the installed
mpmath package. Each command writes to a new external output directory.
The accepted covering check needs only a retained small proposal file:

```text
python research/certify_product_coverings.py --input results/product-symbolic/runs/product-coverings988-scaled2/windows-4.npz --out D:/Codex/new-product-certificate --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages --digits 80
python research/verify_product_coverings.py --source D:/Codex/new-product-certificate --out D:/Codex/new-product-independent --mpmath-path D:/Codex/mpmath-environment/Lib/site-packages --digits 110
python research/audit_product_coverings.py --source D:/Codex/new-product-certificate --out D:/Codex/new-product-primary-audit
python research/audit_product_coverings.py --source D:/Codex/new-product-independent --out D:/Codex/new-product-independent-audit
```

To repeat the numerical proposals, resolve the original Case 988 dense input
from the [numerical campaign](../product-alphabet/README.md), then use:

```text
python research/find_product_cycles.py --source D:/Codex/input/case-988-states.npz --out D:/Codex/new-product-cycles --seeds 32
python research/search_product_returns.py --source D:/Codex/new-product-cycles/cycle-14.npz --out D:/Codex/new-product-returns --radius 0.002 --samples 16385 --steps 18 --max-distance 0.04
python research/refine_product_returns.py --source D:/Codex/new-product-returns --out D:/Codex/new-product-refined --prepend 2 --tails 16
python research/propose_product_coverings.py --source D:/Codex/new-product-refined --out D:/Codex/new-product-windows --per-length 12
```

Cycle ordering and the proposal selection can vary across numerical
environments; the retained accepted windows remove that dependency from
the proof replay. The failed proposals and their exact executed source
versions are preserved in their run directories. Current code contains
the corrected solver scaling and interval expressions.

The figure and instance-binding check use the external dense input:

```text
python research/check_product_symbolic_inputs.py --dense D:/Codex/input/case-988-states.npz --certificate D:/Codex/new-product-certificate --independent D:/Codex/new-product-independent --primary-audit D:/Codex/new-product-primary-audit --independent-audit D:/Codex/new-product-independent-audit --out D:/Codex/new-product-input-check
python research/render_product_coverings.py --input D:/Codex/new-product-certificate/input.json --states D:/Codex/input/case-988-states.npz --metadata results/product-alphabet/runs/product-alphabet-dense01/case-988.json --image docs/assets/diagnostics/new-product-coverings.png --out D:/Codex/new-product-figure
```

`manifest.json` binds every retained packet file, excluding itself. Its scope
does not include later article wording or shared presentation metadata.

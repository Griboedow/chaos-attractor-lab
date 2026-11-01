# Product-posterior orientation certificate

This frozen extension proves a **least-period-92 orbit with negative return
determinant inside the positive-entropy subsystem of Case 988**. It uses the
unchanged [covering packet](../product-symbolic/README.md) and exact parameters
from the [numerical campaign](../product-alphabet/README.md).

The [article, Section 6.1](../../docs/product-alphabet/manuscript.md#61-a-reversed-return-and-an-equivalence-obstruction)
gives the proof. The result excludes ambient-neighborhood topological conjugacy
to orientation-constant local diffeomorphisms, including nonsingular generalized
three-dimensional Henon maps, and smooth time maps and transverse return maps
of Euclidean flows in the corresponding dimensions. Attraction, identification
with the numerical cloud, intrinsic conjugacy, equivalence to known noninvertible
maps and historical uniqueness remain open.

## Accepted certificates

| Run | Result |
| --- | --- |
| `product-orientation988-certificate01` | Factorized posterior at 100 digits; 23 boxes cover all 23 windows. |
| `product-orientation988-independent01` | Full eight-site formula and forward automatic derivatives at 140 digits; 24 boxes. |
| `product-cycles988-certificate03` | Krawczyk existence, least period and real multipliers for cycle identifiers 14, 16 and 23 at 180 digits. |
| `product-cycles988-independent03` | Full-site recomputation at 220 digits and twice the root-box radius. |
| `product-orientation988-audit02` | Exact partitions, binary endpoints, determinant products, Krawczyk arithmetic, root brackets, least-period logic and cross-packet hashes. |
| `product-orientation988-render02` | Reviewed diagnostic schematic and interval plot. |

For G = F^4, determinants are negative on complete windows 4, 5 and 6,
and positive on the other twenty windows. The complete excursion has determinant
in [-3,198,898, -173,698] by the primary formula, independently bounded by
[-5,197,901, -198,607]. All normalization radii exceed 2.3799.
Window 0 and window 1 separation excludes least periods 1, 2 and 4;
the negative return determinant excludes 23 and 46 by the even-power rule.

The short cycles have periods 4, 4 and 7. Only cycle 14 is certified inside
window 0 and can realize the self-loop in the entropy subsystem. No subsystem
or cloud membership is claimed for cycles 16 and 23. The period-92 proof does
not depend on the short-cycle fits.

Both implementations share mpmath interval arithmetic, exact windows and proof
drivers. The full-site mode uses squared distances, eight weights, nested
reflections and automatic derivatives instead of Cartesian factorization and
its analytic derivative. This is formula and derivative independence, not
independent arithmetic software.

## Inventory and unsuccessful work

[Inventory](inventory.json) records every file in 14 runs. [Omissions](omissions.json)
identifies three optional arrays totaling 185,246,700 bytes in the external
hash-addressed cache. Every accepted certificate leaf and exact input is retained.
[Selected summary](selected-summary.json) records the accepted results;
[manifest](manifest.json) binds the frozen packet.

Short-cycle version 1 used a longer exact-input expression; version 2 simplified
it; version 3 exported additional intermediate quantities for the rational audit.
The first audit precedes an explicit cross-packet source-report check. The first
render precedes a legend-position correction. These are retained intermediate
versions, not failed certificates.

The third grid proposal reached 8,451,492 occupied cells and exhausted its
eight-million-new-cell budget. Its corner inflation was sampled, not rigorous.
Two elementary-domain scouts used 1,048,576 and 2,097,152 Sobol states and tested
144 and 80 box/ball domains with normalization cutoffs. All 224 combinations
had numerical violating witnesses. These unsuccessful proposals do not prove
trapping.

[Literature review](literature-review.json) records equations, figures and access
limits. Third-party PDFs and review images remain outside Git, identified by
hash. Four-wing flows, generalized Henon geometries and irregular projected
training are prior art. The addition is the certificate for these exact inputs
and a scoped equivalence obstruction, not a new general orientation theorem.

## Reproduction

Python, NumPy and mpmath 1.3.0 suffice for the accepted interval calculations.
Run these executed source snapshots with fresh external output directories.

```bash
python results/product-orientation/runs/product-orientation988-certificate01/certify_product_orientation.py \
  --source results/product-symbolic/runs/product-coverings988-certificate04 \
  --out /external/orientation-primary --digits 100
python results/product-orientation/runs/product-orientation988-independent01/certify_product_orientation.py \
  --source results/product-symbolic/runs/product-coverings988-certificate04 \
  --out /external/orientation-full --engine full --digits 140
python results/product-orientation/runs/product-cycles988-certificate03/certify_product_cycles.py \
  --source results/product-symbolic/runs/product-cycles988-01 \
  --windows results/product-symbolic/runs/product-coverings988-certificate04/input.json \
  --out /external/cycles-primary --digits 180 --radius 1e-60
python results/product-orientation/runs/product-cycles988-independent03/certify_product_cycles.py \
  --source results/product-symbolic/runs/product-cycles988-01 \
  --windows results/product-symbolic/runs/product-coverings988-certificate04/input.json \
  --primary /external/cycles-primary --out /external/cycles-full \
  --engine full --digits 220 --radius 2e-60
python results/product-orientation/runs/product-orientation988-audit02/audit_product_orientation.py \
  --windows results/product-symbolic/runs/product-coverings988-certificate04/input.json \
  --covering-audit results/product-symbolic/runs/product-coverings988-audit04/report.json \
  --orientation /external/orientation-primary /external/orientation-full \
  --cycles /external/cycles-primary /external/cycles-full \
  --out /external/orientation-audit
```

The covering audit is a dependency, not newly recomputed by the last command.
Its [reproduction instructions](../product-symbolic/README.md) validate the
coverings and full-site recomputation. Scout arrays are unnecessary for these
proofs; source snapshots and seed/range settings permit their regeneration.

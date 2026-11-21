# Toroidal posterior reflections

**Across 3,801,088 distinct parameter/operator configurations, 88 retained forward portraits passed their long numerical checks.** The base survey retained 57. Two targeted generations found 31 more with positive finite-time Lyapunov exponents in all three state directions. All 88 have five million measured binary64 steps, 17 independent starts, a 250,000-point native portrait and an MPFR replay. [External contact-sheet index and checksums](atlas.json) · [All 88 cloud records](manifest.json) · [Three-positive study, two-basin map and 31 cases](three/README.md) · [Exact-input hash census](distinct-configurations.json)

**Sixteen maps have a certified compact positive-entropy subsystem.** The four base examples span all three surface operators, with per-update lower bounds $\log(2)/26$ for Diagonal blades, $\log(2)/104$ for Four curled fins, $\log(2)/295$ for Helical chain and $\log(2)/321$ for Stepped arches. Twelve more certified examples occur among the three-positive cases. Both inverse branches contract, so each certified return contains a uniformly expanding Cantor subsystem conjugate to the full one-sided binary shift. The forward cloud's equality with the corresponding certified subsystem remains open. [Four base proofs](snapback/README.md) · [Twelve follow-up proofs](three/snapback/audit.json)

The base GPU campaign tested 1,572,864 three-state parameter sets. Of 60 visually or dynamically selected long runs, 57 retained positive first-direction growth across ten 500,000-step windows and a majority of 17 starting states. Twenty-three retained second-direction growth above 0.015 in all ten long windows. The 57 are parameter sets. The comparison does not establish 57 different equivalence classes or historical priority. [Case table including exclusions](cases.csv) · [Reproduction manifest](manifest.json)

| Case | Surface operator | Numerical geometry | Five-million-step QR spectrum | Positive starts |
|---:|---|---|---|---:|
| [2189166](case-2189166.json) | Blended linked tubes | Helical chain | 0.260, −0.058, −0.197 | 17/17 |
| [1322937](case-1322937.json) | Nearest linked tube | Slanted frames | 0.336, 0.123, −0.790 | 17/17 |
| [1291273](case-1291273.json) | Nearest linked tube | Layered cusps | 0.539, 0.172, −0.674 | 17/17 |
| [475471](case-475471.json) | Single tube | Forked shell | 0.518, 0.187, −0.637 | 17/17 |
| [1523981](case-1523981.json) | Nearest linked tube | Pleated shield | 0.386, −0.206, −1.844 | 17/17 |
| [370279](case-370279.json) | Single tube | Ribbed fans | 0.158, −0.271, −0.742 | 14/17 |
| [516308](case-516308.json) | Single tube | Domed veil | 0.295, 0.013, −0.931 | 14/17 |
| [2503686](case-2503686.json) | Blended linked tubes | Coiled edge | 0.060, −0.249, −0.615 | 17/17 |
| [2361186](case-2361186.json) | Blended linked tubes | Paired spiral buds | 0.683, −0.351, −1.143 | 17/17 |

“Case” identifies one stored configuration. A positive start has a leading exponent above 0.015 in each of its four windows and passes the finite period test. Different starts can sample different invariant sets. Case 186800 reaches a second rate of 0.377, but its smallest covariance fraction is only 0.00057. The slanted frames and layered cusps retain smallest fractions of 0.202 and 0.194.

## Map and geometry

Let $`b_i\in\mathbb R^3`$ be the equally weighted sites, with temperature $`\tau>0`$. A Gaussian posterior mean replaces one projection of the relaxed alternating-reflection algebra. The earlier repository experiments used an ordinary circle or sphere. This study replaces that surface with one torus, a pair of linked tori, or a smooth blend of their projections. The underlying RAAR algorithm and Gaussian averaging have established antecedents. A targeted literature search did not identify these parameter sets or establish historical uniqueness. [Luke, 2005](https://arxiv.org/abs/math/0405208) · [Previous circle study](../../../docs/soft-constellations.md) · [Previous sphere study](../../../docs/soft-sphere/manuscript.md)

```math
 w_i(x)=\frac{\exp(-\|x-b_i\|^2/\tau)}
               {\sum_j\exp(-\|x-b_j\|^2/\tau)},\qquad
 m(x)=\sum_iw_i(x)b_i,\qquad v=2m(x)-x,
```

```math
 F(x)=\beta x+\beta P(v)+(1-2\beta)m(x),\qquad 0<\beta<1.
```

The first torus has a unit centerline circle in the xy-plane and tube radius $`r`$. Write $`c=(v_x,v_y,0)/\sqrt{v_x^2+v_y^2}`$. Its nearest surface point is $`P_0(v)=c+r(v-c)/\|v-c\|`$ wherever that expression is single-valued. The linked companion has a unit centerline circle in the xz-plane centered at $`(1,0,0)`$. Its projection $`P_1`$ uses the corresponding translated coordinates. The centerline distance is at least one, and the linked runs use $`r<1/2`$, so the tubes are disjoint. For the nearest-surface map, the deterministic tie convention selects the first tube. Derivatives are checked on regular branches.

In the blended map, let $`d_i=\|v-P_i(v)\|`$ and use $`a_i=\exp(-d_i^2/\kappa)/\sum_j\exp(-d_j^2/\kappa)`$. Then $`P(v)=a_0P_0(v)+a_1P_1(v)`$. The blend is smooth across an equal-distance surface when the two individual projections are regular. The numerical and interval routines check the axes and centerlines separately. The exact parameter records store $`1/\tau`$ directly because this coefficient was drawn and rounded in binary32 before the independent replays.

**All forward selections are bounded.** Set $`S=\max_i\|b_i\|`$ and take $`B=1+r`$ for one tube or $`B=2+r`$ for the linked geometries. The same bound holds for mixtures and any choice among tied nearest points. Therefore

```math
 \|F(x)\|\le\beta\|x\|+C,\qquad
 C=\beta B+|1-2\beta|S.
```

Every ball centered at zero with radius strictly above $`C/(1-\beta)`$ absorbs all bounded starting sets. This bound does not identify the limiting invariant subset sampled by a particular portrait.

## Survey and validation

| Surface operator | Screened | Short growth flags | Spatial rescreens | Deep survivors | Shape-eligible | GPU two-rate flags | Independent two-rate survivors |
|---|---:|---:|---:|---:|---:|---:|---:|
| Single ring torus | 524,288 | 41,132 | 25,851 | 23,455 | 18,228 | 113 | 110 |
| Nearest linked tube | 524,288 | 15,333 | 6,154 | 5,700 | 4,658 | 16 | 16 |
| Blended linked tubes | 524,288 | 27,360 | 6,334 | 5,982 | 4,581 | 41 | 41 |

Each of 128 batches per operator contains 4,096 configurations with tetrahedral, octahedral, cuboid, random six-site or perturbed eight-site geometry. All GPU runs used the NVIDIA GeForce GTX 1080. The short stage discards 512 updates and tests four 192-step tangent windows. Spatial covariance screens remove flat or concentrated candidates from the general geometry stream. The deep stage discards 8,192 updates, measures four 1,024-step windows, and records a 768-state terminal cloud. Both linked maps additionally require substantial visits to the region nearer each tube. A 42-coordinate descriptor and farthest-first selection produce 192 inspection cases per operator. [Sampling code](../../../research/shape_search_20261121/projection/screen.py) · [Spatial selection](../../../research/shape_search_20261121/projection/select_shapes.py) · [Checkpoint counts and checksums](checkpoints.json)

A separate search starts from all 35,137 deep survivors. It transports two tangent vectors through four 2,048-step GPU windows and retains cases with second growth above 0.025 in every window. All 170 such flags receive a full native binary64 QR calculation after 32,768 discarded updates and four measured windows of 32,768. The three single-torus exclusions are retained in the audit. [Single-torus audit](two-rate-audit.json) · [Linked-tube audit](linked-two-rate-audit.json) · [Blended audit](smooth-two-rate-audit.json)

The 60 long tests use 100,000 discarded updates, five million measured updates and a full three-vector QR spectrum in ten windows. Each portrait file contains 250,000 consecutive terminal states in its native three coordinates, converted to binary32 for storage after binary64 evolution. Each of 17 starts has a 25,000-step discard followed by four 16,384-step windows. A maximum relative return residual below $`10^{-8}`$ on 768 terminal comparisons flags a period from 1 through 512. Candidates require growth above 0.015 and absence of a detected period in the primary run and at least nine starts, plus a finite-difference Jacobian check. Cases 418637, 1455969 and 1212574 remain in the case table with their actual failed checks.

The [device control](device-control.json) compares short states and tangent growth on the GPU with an independently written CPU map. CPU Jacobians use the two principal tangent directions of the torus and the posterior covariance. Multiprecision replays use three independent MPFR trajectories, reorthogonalizing two displacements every four original updates. Each retained case has the 192-bit run. Additional selected cases have 256-bit runs. A matched 5×5 radius/temperature neighborhood uses factors 0.94, 0.98, 1, 1.02 and 1.06. The saved windows retain regular and competing behavior. [Multiprecision windows](precision/) · [Parameter neighborhoods](neighborhood/) · [Full case and QR table](cases.csv)

## Local periodic dynamics

The base survey stored 202 locally certified periodic orbits from near-return searches and 13 additionally certified attracting companions. Both derivative implementations rechecked all 215 base certificates. A separate 512-bit MPFR calculation checked 1,612 centers across the base and follow-up periodic, sink and entropy packets, using centered finite differences of direct orbits for each inverse Jacobian. The largest center residual is below 8.5×10⁻⁶⁷ and the largest inverse-Jacobian residual row norm is below 5.5×10⁻⁶⁵. [Cross-library numerical audit](mpfr-certificate-audit.json) [Two-evaluator audit](certificate-audit.json)

Each accepted certificate encloses a unique periodic solution, excludes all proper divisors of its period and determines its unstable dimension. The checker uses 110-decimal-place interval arithmetic with automatic differentiation. For the derivative product around the cycle it applies the Cayley transformation to the cubic characteristic polynomial and counts sign changes in the interval-separated first column of its Routh array. Zero or ambiguous pivots fail the certificate. The exact decimal enclosure centers and preconditioners can be fed directly to the verifier. [Cycle certificates](cycles/) · [Attracting companions](sinks/) · [Checker](../../../research/shape_search_20261121/projection/certify_cycles.py)

Cases 1322937 and 1291273 each have two distinct fixed points with unstable dimensions one and two. Case 1322937 has two distinct least-period-two cycles with those same two indices. In the blended Case 2189166, nearby numerical return searches locate cycles with one, two and three unstable directions. The fully repelling period-five orbit supports the certified branching construction. Numerical proximity to a sampled orbit alone cannot establish membership in its attracting set.

Attracting periodic companions explain several observed regular outcomes. Cases 370279 and 516308 have an attracting fixed point. Case 1518001 has a certified attracting orbit of least period 16. Case 175369 has both an attracting fixed point and an attracting period-two orbit. Basin memberships of individual far-away starts are numerical. [Period-16 certificate](sinks/cycles-1518001.json) · [Two attracting periods](sinks/cycles-175369.json)

## Shape comparison

Equal-unit, centered clouds are normalized by their RMS radius. Alpha-complex persistence over the two-element field is measured for the finite offsets of two disjoint temporal halves using samples of 1,536 and 3,072 points per half. Filtration values are converted from squared radii to radii. The comparison takes the maximum of the finite bottleneck distances in dimensions zero, one and two. It retains six intervals per degree and bounds the effect of discarded intervals by their possible matching distance to the diagonal. Unbounded component intervals agree and are omitted. The four matched sample distances are summarized by their median. [GUDHI alpha-complex definition](https://gudhi.inria.fr/python/latest/delaunay_complex_ref.html) · [Controls and per-case diagrams](topology/) · [108-reference comparison](topology-reference.json)

Case 1523981 is the one retained candidate for which the conservative lower bound against every archived reference, 0.0755 RMS units to the displayed precision, exceeds the upper bounds between its own two temporal halves at both sample sizes, 0.0633 and 0.0477. Dadras has the closest truncated diagram median. These finite diagrams describe sampled offsets and cannot distinguish every topological or dynamical equivalence. The general gallery also uses the repository's rotation- and scale-invariant distribution comparisons.

The hard and blended linked maps share the same 524,288 site configurations, tube radii, relaxations and primary starts. Of these pairs, 6,795 received a short growth flag in both maps, 8,538 only with hard selection and 20,565 only with blending. A further 488,390 received neither flag. The [matched long controls](matched-geometry.json) compare growth windows, period outcomes and the repository's shape fingerprints for selected pairs.

## Reproduce

Run from the repository root with Python 3.11. The original environment used NumPy 2.4.6, SciPy 1.17.1, Numba 0.68.1, PyOpenCL 2026.1.5, Pillow 13.0.0, gmpy2 2.3.2, GUDHI 3.13.0 and mpmath 1.4.1. A binary64 portrait replay needs NumPy and Numba. The external cache stores the large arrays. The per-case JSON files in this directory contain the complete sites, parameters and primary starting state.

```powershell
$src = 'research/shape_search_20261121/projection'
$res = 'results/shape_search_20261121/projection'
$cache = 'D:/Codex/chaos-attractor-cache/shape-search-20261121/projection'
python "$src/portrait_atlas.py" --source $res --out "$cache/atlas" --manifest "$cache/atlas-replay.json"
python "$src/reproduce.py" --record "$res/case-2189166.json" --out "$cache/replay-2189166.npz"
python "$src/certify_cycles.py" --verify --source $res --out "$res/cycles" --cases 1322937 1291273 2189166
```

To rebuild the initial survey, run each geometry into its own empty cache subtree. The following sequence reproduces the blended stream. For the other two streams, set $`\mathrm{geometry}`$ to `single` or `linked` and give each command its corresponding cache subtree.

```powershell
python "$src/screen.py" --geometry smooth --cache "$cache/replay-smooth-survey" --batches 128
python "$src/deep.py" --survey "$cache/replay-smooth-survey" --out "$cache/replay-smooth-deep"
python "$src/select_shapes.py" --deep "$cache/replay-smooth-deep" --out "$cache/replay-smooth-select.npz" --limit 192
python "$src/inspect_cases.py" --selection "$cache/replay-smooth-select.npz" --out "$cache/replay-smooth-inspection"
python "$src/second_growth.py" --geometry smooth --deep "$cache/replay-smooth-deep" --out "$cache/replay-smooth-two-growth"
python "$src/audit_two.py" --source "$cache/replay-smooth-two-growth" --cache "$cache/replay-smooth-two-audit" --out "$cache/replay-smooth-two-audit.json"
```

The [manifest](manifest.json) gives the portable cloud path relative to the campaign cache, its SHA-256, the source record, burn-in and measured length. The checkpoint ledger identifies the original cache files for each batch. Long chaotic trajectories depend on floating-point implementations, so cross-device comparisons use finite-time statistics and separate short pointwise controls.

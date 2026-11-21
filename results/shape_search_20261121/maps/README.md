# Three-axis twist maps, 2026-11-21

We continued the repository's [three-axis twist search](../../../docs/twist-geometry.md) with fresh coefficient designs. [The manifest](manifest.json) lists the retained cases, numerical controls and local trajectory paths. Each cloud contains 240,000 consecutive three-state iterates under the key `points`. Iteration used binary64 arithmetic. The saved display coordinates use binary32.

## Retained cases

The per-case links contain parameters and initial states. The manifest gives both an absolute path and a `cloud_relative` path under the shared search cache. Inspect `maps/clouds/<case>.npz` using `numpy.load(path)["points"]`. Phase counts without a certificate describe finite numerical point clouds. Positive starts have positive leading growth in every window.

<!-- case-table-begin -->

| Case and controls | Spatial form | Independent Lyapunov spectrum | Positive starts | Cyclic groups |
| --- | --- | --- | --- | --- |
| [10001029](checks/10001029.json) | Compact folded shell | 0.2060, 0.0285, -0.7475 | 16/16 |  |
| [10003563](checks/10003563.json) | Open oval shell | 0.0898, 0.0305, -0.4234 | 15/16 |  |
| [1001910](checks/1001910.json) | Twisted pennant | 0.2736, -0.1245, -1.2852 | 16/16 |  |
| [10103187](checks/10103187.json) | Bent aperture | 0.0592, 0.0138, -0.3706 | 16/16 |  |
| [10203484](checks/10203484.json) | Parallel pleats | 1.0118, 0.0139, -2.0341 | 16/16 |  |
| [10803912](checks/10803912.json) | Curved mesh canopy | 0.0729, 0.0290, -0.7097 | 16/16 |  |
| [11000087](checks/11000087.json) | Perforated veil | 0.0588, 0.0223, -0.5374 | 16/16 |  |
| [1100117](checks/1100117.json) | Angular frame | 0.3067, -1.0156, -2.2513 | 16/16 |  |
| [1201202](checks/1201202.json) | Vaulted ribbon | 0.1339, -0.1126, -0.4972 | 16/16 |  |
| [1400922](checks/1400922.json) | Parallel hoops | 0.7157, -0.9186, -2.9238 | 16/16 |  |
| [1401913](checks/1401913.json) | Descending arches | 0.3169, -0.2041, -1.0479 | 16/16 |  |
| [1500462](checks/1500462.json) | Folded winglets | 0.1484, -0.3446, -0.6224 | 16/16 |  |
| [1601611](checks/1601611.json) | Scalloped crescent | 0.1856, -0.1467, -0.4898 | 16/16 |  |
| [1700010](checks/1700010.json) | Folded crown | 0.1031, -0.1829, -0.3504 | 16/16 |  |
| [1700146](checks/1700146.json) | Toothed paired arcs | 0.0785, -0.2847, -0.6608 | 16/16 | 2 |
| [1700892](checks/1700892.json) | Woven trough | 0.1567, -0.2053, -0.9242 | 16/16 |  |
| [1701636](checks/1701636.json) | Openwork folded patch | 0.0476, -0.0151, -1.0219 | 16/16 |  |
| [20002293](checks/20002293.json) | Five subgroups in each of 34 phases | 0.0035, -0.2423, -0.7850 | 16/16 | 170 |
| [2003860](checks/2003860.json) | Three narrow return strips | 0.0658, -0.7706, -1.6198 | 16/16 | 3 certified |
| [20100964](checks/20100964.json) | Five subgroups in each of 29 phases | 0.0032, -0.2550, -1.3418 | 12/16 | 145 |
| [20202568](checks/20202568.json) | Six subgroups in each of 13 phases | 0.0074, -0.1481, -0.3613 | 16/16 | 78 |
| [20501481](checks/20501481.json) | Eight folds in each of three phases | 0.0223, -0.1438, -0.3542 | 16/16 | 24 |
| [20600169](checks/20600169.json) | Paired subgroups in 170 phases | 0.0019, -0.2404, -0.7853 | 16/16 | 340 |
| [20601073](checks/20601073.json) | Four-way subdivision of 170 phases | 0.0005, -0.2409, -0.7834 | 16/16 | 680 |
| [20800299](checks/20800299.json) | Seven-way subdivision of 78 phases | 0.0008, -0.1418, -0.3610 | 16/16 | 546 |
| [2100947](checks/2100947.json) | Six separated return threads | 0.0158, -1.1881, -1.7676 | 16/16 | 6 certified |
| [2202255](checks/2202255.json) | Thirty-four-part necklace | 0.0229, -0.2573, -0.7894 | 16/16 | 34 |
| [2303546](checks/2303546.json) | Eight folded strokes | 0.0488, -0.3328, -1.6079 | 10/16 | 8 |
| [2303801](checks/2303801.json) | Four-part return cycle | 0.0934, -0.9361, -1.8173 | 16/16 | 4 certified |
| [2401313](checks/2401313.json) | Five separated return strips | 0.0530, -1.3971, -1.7795 | 16/16 | 5 certified |
| [2503155](checks/2503155.json) | Five-part folded branches | 0.1485, -1.5890, -1.7248 | 16/16 | 5 certified |
| [2601546](checks/2601546.json) | Thirteen-part open ring | 0.0552, -0.1751, -0.3821 | 16/16 | 13 |
| [2602435](checks/2602435.json) | Seven separated return strips | 0.0443, -0.8652, -2.2937 | 16/16 | 7 certified |
| [2603268](checks/2603268.json) | Three folded scrolls | 0.0983, -0.0768, -0.4919 | 16/16 | 3 |
| [2603621](checks/2603621.json) | Twenty-nine-part archipelago | 0.0263, -0.2891, -1.3308 | 12/16 | 29 |
| [2603973](checks/2603973.json) | Paired branching folds | 0.2205, -0.4092, -0.6036 | 16/16 | 2 |
| [2803715](checks/2803715.json) | Three separated return bands | 0.0810, -0.9428, -1.1720 | 16/16 | 3 certified |
| [3201482](checks/3201482.json) | Five curved blades | 0.1053, -0.4673, -1.1854 | 16/16 | 5 |
| [3203938](checks/3203938.json) | Four folded crests | 0.1360, -0.2051, -0.4100 | 16/16 | 4 |
| [3603818](checks/3603818.json) | Eleven-strand comb | 0.0767, -0.2970, -1.3448 | 16/16 | 11 |

<!-- case-table-end -->

## Dynamics and certified cyclic structure

The state first changes to $`r x+(1-r)e_1`$. Three rotations follow, in the ordered planes $`xy`$, $`yz`$ and $`zx`$. In each rotation the current perpendicular coordinate $`w`$ determines the angle $`a+kw+hw^3`$. Each stage fixes $`w`$, so reversing the angle inverts that stage. Reversing the stage order and undoing the contraction gives a global smooth inverse. The Jacobian determinant is $`r^3`$. The [preceding study](../../../docs/twist-priority-and-mechanism.md) identifies an embedded pumped kicked top already studied in 1988. Period doubling and spatially separated chaotic bands have a substantial prior literature. [Pepłowski and Stefański](https://doi.org/10.1016/0375-9601(88)90503-8), [classical damped kicked-top comparisons](https://arxiv.org/abs/2111.13477).

For each successful [trapping certificate](certificates/), the supplied three-dimensional parallelepiped $`T`$ and integer $`m`$ satisfy

```math
F^m(T)\subset\operatorname{int}(T),\qquad F^j(T)\cap T=\varnothing\quad(1\le j<m).
```

Global injectivity implies that the compact sets $`F^i(T)`$ for $`0\le i<m`$ are pairwise disjoint. An intersection for $`i<j`$ would give an intersection of $`T`$ and $`F^{j-i}(T)`$ by applying $`F^{-i}`$. The nested compact connected sets $`F^{mn}(T)`$ have a nonempty connected intersection $`A_0`$. The return map takes $`A_0`$ onto itself. Its images $`A_i=F^i(A_0)`$ are nonempty, compact, connected and pairwise separated. Their union has exactly $`m`$ connected components, and $`F`$ permutes them cyclically. The nested compact sets converge to their intersection in Hausdorff distance. Strict return inclusion therefore gives a neighborhood attracted by this union. Brouwer's fixed-point theorem applied to $`F^m`$ on $`T`$ gives a point of exact period $`m`$ because the earlier images avoid $`T`$.

The proof concerns this localized attracting set. Finite numerical clouds sample its apparent structure. A positive numerical expansion rate alone leaves asymptotic chaos, exact entropy and transitivity unresolved. A rigorous start-entry field, when present, relates the specified initial state to the certified neighborhood.

Each component in a trapping certificate is a cellular continuum in three-space. The sets $`F^{mn}(T)`$ are embedded closed three-balls, and each is contained in the interior of its predecessor. Their intersection is $`A_0`$ by construction. Applying the global diffeomorphism $`F^i`$ gives the same nested-ball description for every $`A_i`$. In particular, each component has trivial positive-degree Čech cohomology. This property concerns the full invariant sets, including fine scales unresolved by a portrait.

The certified attracting sets with three, four, five, six and seven connected components are pairwise nonhomeomorphic across distinct component counts because a homeomorphism preserves that count. This comparison concerns the entire localized invariant sets in the supplied boxes. The [six-component certificate](certificates/2100947-cycle-directed.json) and its [Arb replay](certificates/2100947-cycle-arb.json) locate the canonical start's entry at iterate 124. For Case 2602435, the [seven-component certificate](certificates/2602435-cycle-directed.json) and its [Arb replay](certificates/2602435-cycle-arb.json) include the canonical start's entry at iterate 10.

For Case 2303801, both [Arb](certificates/2303801-parameter-neighborhood.json) and [directed intervals](certificates/2303801-parameter-directed.json) establish the same cyclic enclosure when every one of the ten coefficients varies independently by up to $`2^{-22}`$ about its stored binary64 value. This includes perturbations of the three initially zero cubic coefficients. The original start enters the nominal box on iterate 35. Entry for every parameter tuple in the full block was unresolved.

The primary verifier parses the stored binary64 parameters and box coefficients as exact real numbers, uses Arb interval arithmetic and checks the full dyadic input cube. It encloses the image of every accepted subbox at every intermediate step. The final image lies strictly inside normalized coordinates $`(-1,1)^3`$. For each earlier image at least one coordinate lies strictly outside $`[-1,1]`$. Separating coordinates may vary between subboxes. The independent directed verifier uses outward-rounded binary64 intervals, exact rational coordinate frames and product-rule derivatives. Machin's identity and alternating rational arctangent series bound $`\pi/2`$. Interval Horner evaluation with explicit Taylor remainder bounds encloses sine and cosine. Case 2303801 also passed a [third implementation](certificates/2303801-independent.json) using mpmath interval arithmetic and forward-mode differentiation. [Pointwise controls](certificates/interval-controls.json) compare both implementations with high-precision point values and reject deliberate changes to the return order and twist coefficients.

### Two symbolic branches in Case 2503155

The [directed interval certificate](certificates/2503155-symbolic-directed.json) and a [five-component trapping region](certificates/2503155-cycle-directed.json) give

```math
h_{\mathrm{top}}(F|A)\ge\frac{\log 2}{5}\approx0.1386294
```

per original iterate on the localized attracting set $`A`$. The [affine spatial frame](certificates/2503155-symbolic-box.json) defines two disjoint slabs by their third normalized coordinate. The native widths of the other two coordinates are fixed by that frame. For $`G=F^5`$, the left and right faces of each slab map strictly past opposite ends of the union of the two longitudinal intervals. For every point in either source slab whose image lies within those outer ends, the two transverse output coordinates remain strictly between $`-1`$ and $`1`$. An [exact rational inclusion check](certificates/2503155-containment.json) places both slabs inside the [larger return box](certificates/2503155-box.json).

For a path crossing either source slab, take the portion of its image between the last arrival at the left face and the next arrival at the right face of a chosen target slab, reversing the orientation when needed. Continuity gives a subpath that crosses the target. Transverse containment keeps that subpath inside the target. Iterating this construction realizes every finite binary word. Compactness, applied to longer words with growing past and future lengths, gives bounded full orbits for every doubly infinite itinerary. At return times the full orbit remains inside the larger return box $`T`$. Every such state belongs to each $`F^{5n}(T)`$, hence to $`A_0`$. The two compact slabs have positive separation. Distinct length-$`n`$ words give $`2^n`$ points separated at one of the first $`n`$ return times. This proves the stated entropy bound on $`A`$. It leaves the measure of the symbolic subsystem, transitivity and equality with the numerical cloud open.

The independent numerical spectrum for this parameter tuple is approximately $`(0.14853,-1.58902,-1.72482)`$. All 16 tested starts had positive growth in every window. The numerical orbit visited five phase groups.

The complete localized attracting set has a rigorously noninteger Hausdorff dimension between 1.039 and 1.45, and topological covering dimension one. Its five-component trapping, positive entropy bound, all stated least periods and dimension inequalities hold when all ten coefficients vary independently within $`2^{-24}`$ of their stored values. The canonical start enters the box at iterate six throughout this parameter block. [Dimension proof and explicit enclosures](dimension-proof.md).

The symbolic inequalities passed a [second implementation in Arb](certificates/2503155-symbolic-arb-replay.json). The two nonempty disjoint slabs supply distinct points in $`A_0`$, which is connected and therefore has Hausdorff dimension at least one. Its return homeomorphism has positive entropy. This rules out a finite graph as its underlying topology. A finite-graph homeomorphism permutes the finitely many branch points and endpoints. An iterate preserves every intervening edge and acts monotonically there. A circle homeomorphism also has zero entropy.

**Periodic orbits.** The same face and transverse inequalities imply an orbit of exact period $`5n`$ for every positive integer $`n`$. In particular, there are at least two distinct period-five orbits. Fix a primitive binary word of length $`n`$ and form the product of its $`n`$ assigned slabs. Write each state as $`x_i=(u_i,v_i,w_i)`$ in normalized coordinates and identify indices cyclically. For variable $`u_i`$, use the residual $`G_u(x_i)-u_{i+1}`$ with the sign chosen so that this residual is negative on the lower face and positive on the upper face. The certified strict crossings give these signs for every choice of the other variables. For variable $`v_i`$, use $`v_i-\operatorname{clip}(G_v(x_{i-1}),-1,1)`$, and use the analogous residual for $`w_i`$. Each transverse residual is nonpositive on its lower face and nonnegative on its upper face.

These continuous residuals have a simultaneous zero by the Poincaré–Miranda theorem. Equivalently, subtract the residual vector from the current point, project coordinatewise onto the product rectangle and apply Brouwer's fixed-point theorem. At a fixed point the face signs force each residual to vanish, including at boundary coordinates. The resulting longitudinal equalities put every image within its assigned target interval. The transverse certificate makes the clipping inactive there, giving a genuine periodic orbit of $`G`$. The primitive word forces its least $`G`$-period to equal $`n`$. Its full return orbit lies in $`A_0`$, so the five disjoint cyclic components force the least $`F`$-period to equal $`5n`$. For $`n=1`$ the two constant words give separate period-five cycles because both representatives lie in $`A_0`$ and in disjoint slabs.

Distinct primitive binary words modulo cyclic rotation yield distinct full-map cycles. Thus the number of cycles of least period $`5n`$ is at least the number of primitive binary necklaces of length $`n`$, namely $`n^{-1}\sum_{d\mid n}\mu(d)2^{n/d}`$, where $`\mu`$ is the Möbius function. The [isolated cycles](certificates/2503155-periodic-points.json) additionally have explicit enclosures and multiplier counts. Interval multiple shooting preconditions the cyclic residual equations on the product of small boxes. A contraction mapping of those boxes into their interiors proves existence and uniqueness locally. On the unit circle the quadratic trace term of the three-dimensional characteristic polynomial strictly dominates the sum of the magnitudes of the remaining terms. Rouché's theorem then places two multipliers inside the unit circle and one outside.

## Numerical evidence

All 71 primitive binary necklaces of lengths one through eight have an individually enclosed saddle orbit, with least periods 5 through 40 and exactly one unstable multiplier per orbit. The longer-word existence statement follows from the crossing argument above.

The initial geometry search used four coefficient designs interleaved in each batch. They cover linear twists, cubic-rich twists, general twists and a weaker pair of axes. Each scout discarded 2,048 steps and measured four windows of 1,024 steps. Finite-time growth, occupancy, skewness and covariance informed the diverse review selection. The next search tested phase order through the fraction of spatial variance explained by iterate index modulo $`m`$ for $`2\le m\le36`$. It discarded 3,500 steps and observed 8,192. Independent 80,000-step burn-ins and 160,000-step continuations preceded the visual review. Phase colors, when used, identify consecutive residue classes in the original iteration.

Five parameter designs covered 16,384 general-shape tuples, 81,920 cyclic-order tuples, 49,152 two-exponent tuples, 18,432 perturbations of six cyclic cases and 9,216 perturbations of three of their descendants. All 175,104 full binary64 coefficient vectors are distinct. The manifest builder reconstructs the designs that saved only passing rows, checks each saved vector against its regenerated input index and counts unique tuples across all five designs. It excludes long continuations, repeated initial states and validation perturbations from this total.

The two-exponent screen discarded 3,072 steps and used four windows across the next 8,192. Of 111 settings with second exponent above 0.009 in all four windows, six received full checks. The [multi-start spectra](two-rate-controls/) give 20 windows from 500,000 observed steps at each of 16 starts. Cases 10803912 and 11000087 each had two positive estimates in every window at all sixteen starts. The first local refinement varied one coefficient or several around each of six parent tuples, then tested whether each previously identified phase had a further two- to eight-way split. The 400,000-step continuation retained 100 of 109 inspected proposals. Four received the full long-run controls.

Case 20002293 divides the 34 phases of parent Case 2202255 into five successive subgroups each. All 14,365 pairs of oriented finite-cloud boxes separate, and the smallest separating gap along a tested axis is 0.0394% of the native bounding-box diagonal. The $`170`$-phase labels agree across the held-out tail and eight repeated starts after alignment by a cyclic offset. Its independent leading exponent is approximately $`0.00351`$ per step. Case 20202568 has $`78=13\times6`$ groups with a minimum oriented-box gap of 1.228% of the diagonal. Case 20100964 has 145 groups on the reference orbit and competing outcomes among the other starts. The full spatial point clouds use equal native units. To see the internal dynamics of just one parent phase, take `points[::34]` for Case 20002293 or `points[::13]` for Case 20202568 and label the successive subgroups modulo five or six.

The second local refinement tested 9,216 tuples around the 170-, 145- and 78-group cases. Of 74 inspected proposals, 68 regained their same subgroup count in 1.2-million-step continuations. Cases 20601073, 20800299 and 20600169 have 680 ($`170\times4`$), 546 ($`78\times7`$) and 340 ($`170\times2`$) finite phase groups. Every leading window was positive in the [20-million-step reference trajectories](long-phase-controls/), three-million-step independent QR trajectories and sixteen 1.5-million-step start controls per case. Higher-order 1,015/1,020 proposals had overlapping short-run phase clouds and were left unselected.

The [native binary64 phase controls](phase-binary64/) rebuild the original 240,000 states, confirm their binary32 export exactly and combine them with 500,000 states from each of eight starts using the repository's independent iteration code. Offsets were estimated against the native reference through nearest-neighbor labels. For Cases 20601073 (680) and 20800299 (546), all sampled labels agree after a constant offset at every repeated start. The 4,233,000 and 4,235,868 combined points give positive oriented-box gaps for all 230,860 and 148,785 pairs respectively. The smallest gaps are $`1.54\times10^{-8}`$ and $`3.56\times10^{-6}`$ of the native reference scale. Phase groups at this resolution concern sampled tails. Visualize the 546-group hierarchy at three scales by taking the full trajectory, then every 13th state, then every 78th. Equal axis scaling within each panel preserves each local form.

For each retained shape, the [case packet](checks/) records the binary64 parameters, primary start, 40 disjoint 50,000-step leading-exponent windows after a 100,000-step burn-in, 16 start outcomes and nine nearby or rounded parameter outcomes. Every start discarded 60,000 steps and observed 250,000. Every parameter perturbation discarded 50,000 and observed 200,000. Independently, the original [model and full QR tangent code](../../../research/twist_map.py) discarded 80,000 and observed 600,000 steps. We compared its full spectral sum to $`3\log r`$, and its pointwise analytic tangent to the separately implemented scout. [Selected high-precision checks](precision/) use a scalar implementation and two nearby trajectories rescaled every four steps at 70 and 95 decimal digits.

The [phase checks](phase/) compare a disjoint reference tail and fresh continuations from eight of the 16 initial states to the phase labels in 180,000 training points, allowing a constant cyclic offset for each orbit. A positive gap between all native phase boxes establishes separation for the stored finite cloud at the reported precision. Oriented-box checks additionally test face normals and cross products of the edge directions. The radius and gap measurements use the actual native scale. Negative box gaps permit the boxes to overlap and do not decide point-cloud separation. A combination of coexisting basins and phase-label accuracy below one can indicate a different outcome from another start.

The [geometry table](geometry.json) compares each retained plot to 96 archived twist clouds from the previous gallery. Three PCA-aligned, equally scaled coordinate-plane histograms give a Hellinger distance minimized over axis signs. The same statistic on disjoint pieces of one cloud estimates a sampling baseline. The measure is sensitive to the density of visits and to near-degenerate PCA directions. It provides a finite-catalog comparison. Worldwide visual priority and dynamical equivalence remain open.

An [additional basin survey](additional-basins.json) tested 4,096 initial states at each of eleven fixed parameter tuples, including the original sixteen starts. In the 145-group Case 20100964 and its 29-group parent 2603621, independent interval calculations establish an attracting period-four orbit. The basin scout converged numerically to that alternate cycle from 890 and 898 of the 4,096 starts respectively. Case 10003563 has an interval-certified attracting fixed point alongside the numerically expanding open shell. Every proof supplies its local neighborhood or its cyclic product enclosures. The tested starts and their outcomes have a separate finite observation window.

An [additional audit](all-archive-comparison.json) used all 4,096 archived scouting tuples at equal tail lengths. It compared 3,442 noncollapsed historical clouds and detected 39 phase-separated candidates in the historical finite-time positive-growth screen. The largest recognized order there was 10. These short historical tails and the longer current confirmations have different observation protocols.

Case 2303546 has two distinct interval-certified attracting cycles of least periods two and 24 at exactly the same parameter tuple. Both have open basins, which are disjoint because the limiting cycles differ. For period 24, 192-bit interval multiple shooting isolates one periodic point in each of 24 disjoint phase boxes. On the unit circle the cubic term of the characteristic polynomial strictly dominates the remaining three coefficients, placing all three multipliers inside. The expanding eight-phase reference trajectory is a separate numerical observation at this tuple. [Basin counts and interval certificates](additional-basins.json).

## Reproduction

Run from the repository root using Python with NumPy, SciPy, Numba and Pillow. The primary trapping checks additionally require `python-flint`. The high-precision trajectory controls use `mpmath`. `SHAPE_SEARCH_THREADS` accepts 1 or 2. The default cache is `D:/Codex/chaos-attractor-cache/shape-search-20261121/maps`. The search programs accept `--cache` when a different cache location is needed.

To regenerate an individual portrait cloud directly from its saved tuple and start, run `python research/shape_search_20261121/maps/replay.py 2503155`. The output goes to the cache's `replay` directory and the program reports the new and archived hashes. The complete search commands follow.

```powershell
$code = "research/shape_search_20261121/maps"
$results = "results/shape_search_20261121/maps"
$cache = "D:/Codex/chaos-attractor-cache/shape-search-20261121/maps"
$env:NUMBA_CACHE_DIR = "$cache/jit"
python "$code/search.py" --count 2048 --batches 8
python "$code/preview.py" --batches 0 1 2 3 4 5 6 7 --limit 288 --name preview-initial
python "$code/validate.py" --preview "$cache/preview-initial/preview.json" --proposals "$results/proposals-initial.json"
python "$code/cyclic.py" --first 8 --batches 20 --count 4096
python "$code/phase_review.py"
python "$code/validate.py" --preview "$cache/review-cyclic/preview.json" --proposals "$results/proposals-cyclic.json"
python "$code/validate.py" --preview "$cache/review-cyclic/preview.json" --proposals "$results/proposals-certified.json"
python "$code/validate.py" --preview "$cache/review-cyclic/preview.json" --proposals "$results/proposals-symbolic.json"
python "$code/hyper.py" --batches 12 --count 4096
python "$code/preview.py" --cache "$cache/hyper" --batches 0 1 2 3 4 5 6 7 8 9 10 11 --name review --limit 192
python "$code/validate.py" --preview "$cache/hyper/review/preview.json" --proposals "$results/proposals-hyper.json"
python "$code/two_rate_controls.py"
python "$code/phase_refine.py" --count 3072
python "$code/validate.py" --preview "$cache/nested-phase/preview.json" --proposals "$results/proposals-nested.json"
python "$code/phase_refine.py" --groups 6 7 8 --name nested-deep --count 3072
python "$code/validate.py" --preview "$cache/nested-deep/preview.json" --proposals "$results/proposals-deep.json"
python "$code/phase_validate.py"
python "$code/phase_validate.py" --orders "$cache/nested-phase/preview.json" --cases 20002293 20100964 20202568 20501481
python "$code/phase_validate.py" --orders "$cache/nested-deep/preview.json" --cases 20601073 20800299 20600169
python "$code/long_phase_controls.py" --cases 20601073 20800299 20600169
python "$code/phase_precision.py" --cases 20601073 20800299 20600169 20002293
python "$code/geometry.py"
python "$code/archive_audit.py"
python "$code/precision.py" --cases 1700146 2603268 2202255 1701636 3203938 1400922 2601546 3603818 3201482 1601611
python "$code/precision.py" --cases 20002293 20100964 --observed 64000 --burn 16000
python "$code/precision.py" --cases 20601073 20800299 --observed 128000 --burn 32000
python "$code/basin_scout.py" --count 4096 --cases 2503155 2303801 2602435 2401313 2803715 2100947 2003860 20100964 2603621 2303546 10003563
python "$code/sink_cycles.py"
python "$code/summarize.py"
python "$code/audit_packet.py"
```

The geometry command regenerates the historical parameter grid when the archived dense packet is unavailable. Pass `--archive <dense/data.npz>` to use the original file. Its SHA-256 is checked. The original 96-gallery packet was used for the reported comparison.

To replay one rigorous certificate, run these commands. The second verifier checks the entire stored tree and independently subdivides any unresolved leaf.

```powershell
$cert = "$results/certificates"
python "$code/verify_return.py" "$cert/2303801-box.json"
python "$code/directed_return.py" "$cert/2303801-box.json" "$cert/2303801-verification.json"
python "$code/independent_return.py" "$cert/2303801-box.json" "$cert/2303801-verification.json" --extra 50000
python "$code/verify_return.py" "$cert/2303801-box.json" --parameter-radius-bits 22
python "$code/directed_return.py" "$cert/2303801-box.json" "$cert/2303801-parameter-neighborhood.json" --parameter-radius-bits 22
python "$code/directed_cover.py" "$cert/2503155-symbolic-box.json" "$cert/2503155-faces.json" --output "$cache/2503155-symbolic-repeat.json"
python "$code/arb_replay.py" "$cert/2503155-symbolic-box.json" "$cert/2503155-symbolic-directed.json" --output "$cache/2503155-symbolic-arb-repeat.json"
python "$code/directed_cycle.py" "$cert/2503155-box.json" --output "$cache/2503155-trap-repeat.json" --budget 1250000
python "$code/arb_replay.py" "$cert/2503155-box.json" "$cert/2503155-cycle-directed.json" --output "$cache/2503155-trap-arb-repeat.json"
python "$code/containment.py" "$cert/2503155-symbolic-box.json" "$cert/2503155-box.json" "$cert/2503155-symbolic-directed.json" --output "$cache/2503155-box-inclusion-repeat.json"
python "$code/periodic_points.py" --maxword 8
python "$code/directed_cover.py" "$cert/2503155-symbolic-box.json" "$cert/2503155-faces.json" --primary "$cert/2503155-symbolic-directed.json" --parameter-radius-bits 24 --output "$cache/2503155-uniform-cover.json"
python "$code/directed_cycle.py" "$cert/2503155-box.json" --parameter-radius-bits 24 --budget 300000 --output "$cache/2503155-uniform-return.json"
python "$code/arb_replay.py" "$cert/2503155-symbolic-box.json" "$cert/2503155-symbolic-parameter.json" --output "$cache/2503155-uniform-cover-arb.json"
python "$code/arb_replay.py" "$cert/2503155-box.json" "$cert/2503155-cycle-parameter.json" --output "$cache/2503155-uniform-return-arb.json"
python "$code/slice_bound.py" "$cert/2503155-box.json" "$cert/2503155-symbolic-parameter-arb.json" "$cert/2503155-cycle-parameter-arb.json" "$cert/2503155-symbolic-box.json" --output "$cache/2503155-global-lower.json"
python "$code/slice_geometry.py" "$cert/2503155-box.json" "$cert/2503155-symbolic-parameter-arb.json" "$cert/2503155-cycle-parameter-arb.json" --grid-power 6 --output "$cache/2503155-phase-lower-arb.json"
python "$code/directed_slices.py" "$cert/2503155-box.json" "$cert/2503155-cycle-parameter-arb.json" --grid-power 6 --output "$cache/2503155-phase-lower-directed.json"
python "$code/dimension.py" "$cert/2503155-box.json" --fraction 9 20 --parameter-radius-bits 22 --output "$cache/2503155-upper-arb.json"
python "$code/directed_dimension.py" "$cert/2503155-box.json" "$cert/2503155-dimension-arb.json" --budget 600000 --output "$cache/2503155-upper-directed.json"
```

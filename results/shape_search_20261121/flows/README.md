# Three-dimensional continuous flows

Each portrait uses the three native coordinates of one autonomous ODE. The parameter search used deterministic integer case IDs, double-precision integration and a seed of 2026112123. The initial survey required a largest finite-time variational growth rate above 0.008 in every observed window and a smallest-to-largest spatial covariance eigenvalue ratio above 0.008. The [screen accounting](accounting.json) records each completed batch, including failed cases. A selected parameter set is a numerical chaos candidate. Similar-looking forms occur in the literature, and the finite survey does not settle historical priority.

[Case parameters, cloud hashes and absolute local paths](manifest.json) · [Initial-condition and adaptive-solver controls](controls.json) · [Occupation and angular measurements](geometry.json) · [Diagonal equilibrium enclosures](diagonal-equilibria.json) · [Run settings](run-settings.json)

The search evaluated **70,356 distinct coefficient vectors** after counting repeated initial conditions, step sizes and solvers once. It retained 35 numerically chaotic clouds and 174 numerically attracting periodic knots with 31 different normalized Alexander polynomials. Four collections contain linked coexisting orbits. Two chaotic cases have independently tracked unstable knotted cores. The comparison covers all 215 geometries against 177 archived samples and separately against the 98 samples identified as ODE references. [Unique-equation and control counts](index-summary.json) · [174 portable cloud records](periodic-clouds.jsonl) · [Full archive comparison](archive-comparison.json.gz) · [ODE-only comparison](ode-reference-comparison.json.gz)

| Family | Unique screened equations | Selected chaotic clouds | Validated stable knots |
| --- | ---: | ---: | ---: |
| Cartesian harmonic coupling | 10,560 | 10 | 124 |
| Angular resonance | 13,716 | 10 | 37 |
| Quartic precession | 20,480 | 8 | 13 |
| Smooth switching planes | 24,576 | 3 | 0 |
| Crossed angular resonance | 1,024 | 4 | 0 |

The periodic-knot screen covered Cartesian, angular and precession flows. The 70,356 count excludes the diagnostic bisection probes used to locate a multiplier crossing. Dynamically conjugate equations can have different coefficient vectors. Each archive comparison uses both similarity-normalized and affine-whitened shape descriptors. In Case 180844 the nearest of the 98 ODE references was Lorenz84 in both modes, at distances 0.415 and 0.176 versus 0.023 and 0.023 for two disjoint temporal halves of the chaotic candidate. A finite atlas and these geometry statistics cannot establish priority. [Archive inputs](archive-inputs.json) · [ODE-only rankings](ode-reference-comparison.json.gz)

## Knotted bands and periodic orbits

**Case 180844.** A stable, seven-crossing sampled knot in Case 28756 loses transverse stability as the axial drive coefficient $`e`$ increases from −0.4297716348 to −0.3538341348. At the endpoint, both the numerically chaotic trajectory and the unstable periodic orbit solve the same autonomous equations at identical parameter values. Shooting gives period 26.6537484684, return error $`2.86\times10^{-12}`$ in native coordinates and unstable multiplier −47.00369811. The independent three-exponent estimate is (0.09377, 0.00041, −1.95504). All five tested starts passed the growth screen. Two longer traces retained positive growth in all 40 windows of 5,000 time units after individual transients of 1,200 units. The window rates range from 0.08797 to 0.09466. [Branch and shooting](band-continuation.json) · [Solver, geometry and lifetime controls](band-audit.json)

The fixed cloud comprises 266,666 native XYZ points with RMS radius 0.891452309 about its own arithmetic mean. Distances use this common radius, with no registration or rescaling of individual coordinates. For the original sample, the 95th percentile distance from every fourth skeleton vertex to its nearest cloud point is 0.1122% RMS. The reverse measure from every sixteenth cloud point to the nearest of 32,768 skeleton vertices is 3.8608% RMS. An independent subsampling and exact closest-segment calculation gave 3.8970%. Two fresh state-only trajectories, one DOP853 and one RK45, gave 3.8924% and 3.8787%. Each visited all 64 equal-time bins of the skeleton. The largest reverse distance in these subsamples was about 8% RMS. These finite cloud comparisons describe a band around the periodic trajectory. They do not identify the topology or invariance of the entire chaotic limit set. [Geometric audit](band-audit.json) · [Point and orbit records](manifest.json)

Two independent centered differences of the physical section return near one full period converged to unstable multipliers −47.0036976 and −47.0036934. At a displacement of $`3\times10^{-8}`$ skeleton RMS the Euclidean condition number of $`DP-I`$ was 289.9953. The central return residual was below $`1.34\times10^{-10}`$ skeleton RMS. The near-zero contracting multiplier is poorly conditioned in a computed two-by-two spectrum. The separately integrated analytic divergence gives a determinant of $`2.97\times10^{-22}`$ for the return over one period. [Section refinement, matrices and tolerances](band-audit.json)

**Case 200949** supplies a second parameter path from the stable knot of Case 28066. The continuation gives a real unstable multiplier near −29.15 and a 95th-percentile cloud-to-skeleton distance of 7.86% RMS. Five starts passed its growth screen. The sampled core and precursor have the same normalized Alexander polynomial. These two paths are examples within one angular-resonance family. [Case records](manifest.json) · [Continuation](band-continuation.json)

### A sequence of period doublings

The primary periodic orbit's nontrivial real multiplier crosses −1 near $`e=-0.37574761`$. Adaptive shooting and an independent state return detect attracting periods approximately 1, 2, 4, 8 and 16 times the primary period at the following five points. The quoted bracket in the attached file records numerical bisection. A validated ODE interval enclosure of the flip remains open. [Parameter sweep](period-doubling.json) · [16-period return](long-doubling.json)

| Axial coefficient $`e`$ | Period relative to primary | Refined period | Dominant transverse multiplier |
| ---: | ---: | ---: | ---: |
| −0.376615385 | 1.0000 | 26.62991 | −0.82207 |
| −0.375096635 | 2.0001 | 53.26562 | 0.54397 |
| −0.360668510 | 4.0045 | 106.70534 | 0.02264 |
| −0.358010697 | 8.0098 | 213.45412 | 0.00000909 |
| −0.357631010 | 16.0198 | 426.92341 | 0.79069 |

At $`e=-0.375096635`$, the primary orbit is already unstable with multiplier −1.15481. Its attracting period-two daughter winds twice around its nearest phase, with strictly positive sampled phase increments. Exact polygonal linking has magnitude 21. The daughter's computed Alexander polynomial is the product of the $`T(2,21)`$ polynomial and the primary polynomial evaluated at $`t^2`$. These independent results support a $`(2,21)`$ cable interpretation. [Daughter polynomial](doubling-witnesses) · [Local winding and product check](doubling-cable.json) · [Pair clouds and multipliers](doubling-pairs.json)

At the next two bifurcations, the exact polygonal linking magnitudes of the period-two and period-four pair and the period-four and period-eight pair are 91 and 357. The high-period curves pass extremely close together. At 640 vertices the second pair gave 170 and the intermediate 1,280 mesh gave −275 for the chosen orientation. The signed result stabilized at −357 across the 5,120, 10,240, 20,480 and 32,768 arclength meshes and the original 32,768-point time mesh. The distinct camera projections at a fixed resolution agree even when the mesh is too coarse to resolve the physical strands. [Signed multi-resolution results](doubling-link-refinement.json)

The four observed doubling pairs have absolute linking numbers 21, 91, 357 and 1,435. For the 8-to-16 pair at $`e=-0.357631010`$, the parent multiplier is −1.05912 and the daughter multiplier is +0.79069. The daughter's 95th-percentile distance to the parent is 0.0687% of the parent RMS radius. Three independent initial states produced the long near-return. The fourth pair required much finer resolution: the 32,768-vertex arclength polygon gave −1,439, while the 49,152-vertex polygon and both 65,536-vertex time and arclength polygons gave −1,435 in every tested camera. [Pair clouds](doubling-pairs.json) · [All meshes](doubling-link-refinement.json) · [Exact binary64 pair links](doubling-exact-links.json) · [Period-sixteen return](long-doubling.json)

### Five unstable periodic orbits at one chaotic parameter

Continuing the daughter branches to the final Case 180844 coefficients gives four additional unstable periodic orbits. For order 2, direct full-period shooting gave period 53.3884678 and multiplier +220.08664. For higher orders, the $`C_7`$ symmetry makes each orbit a rotated copy of itself after one seventh of its time period. A relative shooting solve closes the endpoint onto the rotated starting point. Seven rotated sectors reconstruct the trajectory. The full transverse multiplier is the seventh power of the dominant multiplier of the return after undoing the rotation. Two independent state-only adaptive solvers checked each sector closure. [Full portable core records](nested-skeletons.json) · [Relative-shooting checks and continuation](relative-core-validation.json) · [Exact pairwise linking](nested-core-links.json)

| Period order | Period at Case 180844 | Dominant full-period multiplier | Cloud-to-orbit 95th percentile, RMS |
| ---: | ---: | ---: | ---: |
| 1 | 26.65375 | −47.00370 | 3.86% |
| 2 | 53.38847 | +220.08664 | 2.35% |
| 4 | 106.75931 | −39,957.94 | 1.09% |
| 8 | 213.53786 | $`-1.971\times10^9`$ | 0.855% |
| 16 | 427.05819 | $`-3.206\times10^{18}`$ | 0.827% |

All five rows use the same equation and the same chaotic reference cloud. The largest multiplier arises from seven replications of a sector multiplier of −440.26087. The order-16 rotated endpoint discrepancies are $`1.47\times10^{-9}`$ and $`1.16\times10^{-10}`$ in native units. Nineteen-phase measurements test divisors from 2 through 256. The order-16 orbit had a minimum maximum-phase separation of 0.0570 of its own orbit RMS under those fraction-of-period shifts. [Orbit records and scale conventions](nested-skeletons.json)

At the final chaotic coefficients, the exact stored time polygons give signed links −21, −90, −357 and −1,435 between successive orders. The 2T/4T result differs by one from the earlier parameter-pair result. The same signed value of −90 persisted at 65,536 and 131,072 vertices per orbit in two integer projections. A numerical continuation and a discrete linking number alone do not establish topological preservation along the entire smooth branch. [Final-pair exact predicates](nested-core-links.json) · [Final 2T/4T mesh refinement](nested-link-refinement.json)

### Recognition of the full sampled knots

An exact elementary move replaces two polygon edges by the opposite side of their triangle after testing that its interior and boundary meet no other edge except at permitted endpoints. Exact spatial bounding boxes, integer orientation predicates and a second check of the final polygon support the reduction. Cases 180844 and 28756 reduced from 32,768 vertices to short polygons with seven-crossing diagrams. For a further classification, a stored integer axis projection puts every oriented edge at strictly positive angular increment and gives winding two. A closed two-strand braid is a $`T(2,q)`$ torus link, and its determinant is $`|q|`$ for odd $`q`$. The reductions and axis witnesses identify both cores (180844 and 200949) and the stable precursor 28756 as $`T(2,7)`$ up to mirror reflection. The same check gives $`T(2,5)`$ for 21038, and $`T(2,3)`$ for 27762 and for the representative of 26471. The assertion concerns the sampled polygons. The accuracy of their correspondence to the smooth ODE orbit remains numerical. [Dense-polygon reductions](dense-polygon-isotopies.json.gz) · [Two-strand witnesses](two-strand-witnesses.json) · [Kassel and Turaev, *Braid Groups*, section 2.3](https://web.math.ucsb.edu/~bigelow/books/kasselturaev.pdf)

The complete audit used the exact binary64 values stored in the finite XYZ samples. All coordinates of one polygon were multiplied by their common power-of-two denominator before integer reduction. Every one of the 174 stable polygons, both unstable cores and the 2T daughter passed (177 of 177). An independent Kauffman-bracket calculation gives the Jones polynomial in two projections of each of the 174 stable polygons and passes known trefoil, cinquefoil and figure-eight diagrams. The results distinguish 31 signatures after identifying mirror images and 44 handed Jones polynomials. The coefficients in each stored term are paired with actual integer powers of $`t`$, with unknot normalization one. Numerically simulated trajectories joined by line segments and exact algebraic invariants of those polylines represent distinct levels of evidence. [Dense-polygon audits](dense-polygon-isotopies.json.gz) · [Jones polynomials and projection controls](jones-polynomials.json.gz)

### Coexisting linked attractors

A single equation and fixed coefficients give each row below. Equivariance constructs the other cycles from one numerically attracting solution. The shapes have finite separations; a trajectory from one initial condition visits one attractor after settling. The linking number is computed between distinct oriented components, with its absolute value reported below.

| Case | Components | Absolute linking number of every pair | Trials captured by any member |
| --- | ---: | ---: | ---: |
| 2479 | 3 | 30 | 15/15 |
| 26471 | 4 | 6 | 18/18 |
| 23599 | 2 | 10 | 12/12 |
| 25427 | 2 | 12 | 12/12 |

At Case 2479 the dominant transverse multiplier magnitude from independently perturbed section maps was 0.05454–0.05465. At Case 26471 it was $`1.4522\times10^{-6}`$–$`1.4536\times10^{-6}`$. Nearby starts sometimes settled on another symmetry copy in Case 2479. The 640- and 1,280-vertex polygons agreed in three integer projections per resolution. The repository stores the smaller polygons, exact integer crossing counts, a hash of the finer witness and the per-start assignment. A Hopf-link control gives absolute linking number one and a separated pair gives zero. [Link witnesses](link-witnesses) · [Basin assignments and independent Floquet checks](link-checks.json)

### Stable-knot screen

Stable near-returns were refined by shooting in a coordinate section, then independently integrated with DOP853. Each accepted cycle closed within $`5\times10^{-7}`$ RMS and had a numerically neutral autonomous multiplier and both transverse moduli below 0.995. Perturbation checks used 1% and 4% RMS displacements. Comparison of the orbit at 19 phases with shifts of $`T/k`$, for every integer $`k=2,\ldots,256`$, tested repeated traversal in the accepted catalog.

Twelve camera directions supplied diagrams at each of 320, 640 and 1,280 arclength vertices. Cases with an inconsistent determinant across usable projections, or inconsistent normalized Alexander coefficients in the two lowest-crossing usable views per mesh, were excluded. Different polynomials distinguish sampled knots; matching polynomials can belong to different knot types. For the curated examples, projected intersections, orientation signs, depth order and polynomial determinants were recomputed with exact arithmetic on common-RMS integer polygons. Circle, trefoil, cinquefoil and figure-eight controls accompany the algorithm. The exact statements concern these fixed polygons. ODE existence, stability and the isotopy from the smooth orbit to its polygon are inferred numerically. Prior constructions include stable knots and links in ODEs. [All validated cases](periodic-catalog.csv) · [Screen counts and distinct polynomials](periodic-accounting.json) · [Curated polygon witnesses](knot-witnesses) · [Earlier stable-knot construction](https://applmath.cjoe.ac.cn/jweb_yysxxb/EN/10.12387/C1991001)

Cases 27762, 21038 and 28756 have diagrams with three, five and seven crossings. Case 967 has determinant 243 and polynomial coefficients (9, −30, 52, −61, 52, −30, 9). Cases 973 and 2160 both have determinant 49 but different polynomials. Projection crossing counts in the catalog are upper bounds on a knot's minimum crossing number. [Local controls](knot-controls.json) · [Shooting seeds](knot-starts.jsonl)

## Other chaotic forms

Cases 132, 306 and 108 show a six-loop rosette, twin towers and paired diamond folds in the native Cartesian flow. All three passed four broad starts and sixteen local reinjections. Cases 21253 and 20897 have four and eight folded eyes in the angular-resonance flow. The serrated crest in precession Case 52506 passed one of five broad starts, all six local starts and sixteen late perturbations. In three Cartesian examples, interval arithmetic encloses coexisting stable equilibria on the invariant diagonal. The packet preserves the less successful starts and long-screen failures alongside the selections. [Local reinjections](recovery.json) · [Periodic returns within chaotic clouds](periodic-orbits.json) · [Long-screen outcomes](long-screen.csv) · [Interval equilibria](diagonal-equilibria.json)

## Equations

### Cartesian harmonic coupling

For cyclic indices $`(i,j,k)`$, each row of the parameter matrix contains $`(d,a,b,c,\phi,\psi,\chi,m,n)`$ in that order. The components obey

```math
\dot x_i=-d_i x_i+a_i\sin(x_j+\phi_i)+b_i\sin(x_k+\psi_i)+c_i\sin(m_i x_j+n_i x_k+\chi_i).
```

The divergence is the constant $`-\sum_i d_i`$. Each coordinate has the limiting bound $`(|a_i|+|b_i|+|c_i|)/d_i`$. Equal parameter rows give a threefold cyclic symmetry, and zero phases give central inversion. The Thomas cyclic sine flow occurs when $`b=c=0`$ with equal damping and unit forward coupling. Three-dimensional labyrinth chaos is well documented. [Thomas, 1999](https://doi.org/10.1142/S0218127499001383) · [Sprott, Labyrinth chaos](https://sprott.physics.wisc.edu/pubs/paper302.pdf)

The diagonal $`x=y=z=s`$ is invariant when the parameter rows are equal. At a scalar zero on this line, let $`G_j,G_k`$ be the two cross derivatives of one component. The real eigenvalue is $`-d+G_j+G_k`$. The other two have real part $`-d-(G_j+G_k)/2`$ and imaginary parts of opposite sign with magnitude $`\sqrt{3}|G_j-G_k|/2`$. The interval calculation establishes signs throughout each reported root enclosure. It uses the exact binary64 parameter values, checks opposite endpoint signs and bounds the scalar derivative away from zero.

### Angular resonance

Write $`u=x+iy`$, $`s=|u|^2`$, $`h=2/(1+s)`$. The eleven parameters, in archive order, are $`n,\alpha,\beta,\omega,\kappa,\tau,A,d,e,f,\phi`$.

```math
\begin{aligned}
\dot u &= u[\alpha-s+\beta z+i(\omega+\kappa z+\tau s)]
          +A h^{(n-2)/2}\overline{u}^{n-1},\\
\dot z &= -d z+e\frac{s-1}{1+s}
          + f h^{n/2}\operatorname{Re}(e^{i\phi}u^n).
\end{aligned}
```

Rotation by $`2\pi/n`$ about the $`z`$ axis preserves the equations. The axis contains an equilibrium at $`z=-e/d`$. The $`z`$ drive has magnitude at most $`|e|+2^{n/2}|f|`$, giving an eventual axial bound. The $`-s u`$ term dominates planar radial growth at large radius once $`z`$ is bounded. The number of rotational sectors belongs to the field design. Attracting folds, local returns and occupation across sectors depend on integration and the selected parameters. Angular statistics in the result packet use the stored Cartesian trace. Earlier multiwing constructions also use rotational symmetry to regulate wing count. [Yan et al., tunable multiwing construction](https://doi.org/10.1016/j.aej.2024.09.011)

### Forced quartic precession

Let $`v=(x,y,z)`$ and let $`D`$ be positive diagonal. The sixteen coefficients comprise three damping entries, three components of the constant torque $`t`$, then $`a_i,q_i,h_i`$ and $`c`$.

```math
H(v)=a\cdot v+\tfrac12\sum_i q_i v_i^2+\tfrac14\sum_i h_i v_i^4+cxyz,\qquad
\dot v=t-Dv+\frac{v\times\nabla H(v)}{1+|v|^2}.
```

The angular term is orthogonal to $`v`$ and has zero divergence. The divergence of the complete system is $`-\operatorname{tr}D`$, and the radius has limiting upper bound $`|t|/\min_i D_{ii}`$. This model is a constructed mathematical precession flow. Forced dissipative gyrostats and their chaotic solutions have substantial prior work. [Seshadri and Lakshmivarahan, 2023](https://arxiv.org/abs/2307.12337)

### Smooth switching planes

With hidden rows $`(a_r,b_r,w_r)`$ containing three normal components, one offset and three output components, the state obeys

```math
\dot x=t-Dx+\sum_r w_r\tanh(a_r\cdot x+b_r).
```

The flat parameter array begins with the row count, the three damping rates and the three drift components, followed by the seven entries of each hidden row. The bounded switches and positive damping bound every coordinate. The generator fits the local linearization at the origin to an index-one or index-two saddle-focus. These flows are continuous-time recurrent networks. A fitted local saddle-focus alone gives no conclusion about chaos.

### Crossed angular resonance

Two angular fields share the same three state variables. The second has an independently sampled axis and center. If $`F_a,F_b`$ are the angular fields, $`R`$ is the orthogonal rotation, $`c`$ the displaced center and $`w`$ the mixing fraction, the combined velocity is

```math
\dot v=(1-w)F_a(v)+w R F_b(R^{\mathsf T}(v-c)).
```

The archived flat array stores the two 11-parameter vectors, the nine row-major entries of $`R`$, the three components of $`c`$ and $`w`$. Cases 120055, 120102, 120824 and 120424 survived adaptive and initial-condition checks with their recorded restrictions. This variation reuses the single-axis resonance mechanism.

## Numerical evidence

The fixed-step solver integrates the analytic tangent equation together with the state by classical fourth-order Runge–Kutta. After transient removal it renormalizes the tangent every five steps and at window boundaries. The dense portrait files contain consecutive subsampled states; their parameter arrays and SHA-256 hashes are registered in the manifest.

Adaptive replay uses DOP853, a separately coded vector field and centered-difference Jacobians. Integrating the full fundamental matrix with QR normalization estimates all three exponents. A thirteenth state accumulates the divergence for a check against the sum of the exponents. Each control lists the integrated time, solver tolerances, growth windows and checks against the analytic implementation.

To compare occupation across starts the same 32 fixed one-dimensional projections, coordinates, anchor center and scale are used for every run of a case. Geometry distances are mean absolute quantile differences on the 1st through 99th percentiles, expressed in anchor RMS radii. The distance between the two temporal halves of the anchor supplies a sampling comparison. A positive growth rate at another start by itself does not identify the same occupied set.

## Reproduction

CPython 3.12.14, NumPy 2.5.4, SciPy 1.18.1 and Numba 0.68.0 produced the reference arrays. Pillow supports contact sheets. The equilibrium enclosure uses mpmath 1.3.0. The polygon checks support optional python-flint 0.9.0 for integer determinants. From the repository root, choose fresh output files and run

```powershell
$flow = 'research/shape_search_20261121/flows'
$packet = 'results/shape_search_20261121/flows'
python "$flow/replay_case.py" --case 132 --out D:/attractor-replay/flow132.npz
python "$flow/replay_case.py" --case 180844 --out D:/attractor-replay/chaotic-180844.npz
python "$flow/replay_case.py" --case 180844 --skeleton --out D:/attractor-replay/orbit-180844.npz
python "$flow/replay_case.py" --case 26471 --linked --out D:/attractor-replay/linked-26471.npz
python "$flow/replay_case.py" --case 27762 --out D:/attractor-replay/stable-27762.npz
python "$flow/replay_case.py" --case 180844 --skeleton-order 16 --out D:/attractor-replay/final-core-16.npz
python "$flow/replay_case.py" --case 180844 --doubling 16 --out D:/attractor-replay/doubling-pair-16.npz
python "$flow/braid_axis.py" --verify "$packet/two-strand-witnesses.json" --reductions "$packet/dense-polygon-isotopies.json.gz"
python "$flow/knot_witness.py" --verify "$packet/knot-witnesses"
python "$flow/knot_witness.py" --verify "$packet/band-witnesses"
python "$flow/symmetry_links.py" --verify "$packet/link-witnesses"
python "$flow/diagonal.py" --manifest "$packet/manifest.json" --out D:/attractor-replay/diagonal-equilibria.json
```

The replay accepts any case in the manifest or stable-knot catalog. Its printed SHA-256 can be compared with each archived reference. Bitwise equality depends on the numerical runtime. Cloud paths in the manifest are both absolute and relative to the campaign cache root. Every XYZ archive exposes a float64 array named `points`. Linked sets additionally expose `component` and `offsets` so that polylines belonging to different orbits can be drawn separately. The burst geometry and exponent windows vary at individual points across chaotic replays. [Reproduction code](../../../research/shape_search_20261121/flows) · [Run settings](run-settings.json)

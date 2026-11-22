# Cyclic polynomial maps with nested trapping annuli

A single orbit draws each portrait. The four certified maps have local attractors with exactly two, three, four or five connected components. Each component winds around the origin in the sense that it separates the origin from infinity. The three- and five-component permutations also force positive topological entropy on compact invariant sets in the radial gaps. The latter theorem concerns trajectories between the attracting bands. The Lyapunov rates measured in the bands are numerical.

| Case | Rotation order | Certified radial cycle, innermost = A | Mean numerical rates per iterate | Portrait |
| --- | ---: | --- | --- | --- |
| 13745 | 17 | A → B → A | +0.125195, +0.005589 | [Two nested crowns](cyclic-13745.png) |
| 205188 | 18 | A → B → C → A | +0.086178, +0.000386 | [Three nested rings](cyclic-205188.png) |
| 5529 | 18 | A → C → B → D → A | +0.014771, −0.001591 | [Four interleaved rings](cyclic-5529.png) |
| 202080 | 20 | A → C → D → B → E → A | +0.027593, −0.000218 | [Five nested rings](cyclic-202080.png) |

**Drawing date:** 2026-11-22 for the seven portraits in this packet. The other three selected settings are numerical. [Case 1723](cyclic-1723.png) has a small four-cornered central web surrounded by four arcing outer branches. [Case 14725](cyclic-14725.png) has eight concave sides and a large hollow center. [Case 11450](cyclic-11450.png) traces three bent chains of small loops. The second rate in Case 205188 crosses zero in at least one measured window.

## The map and parameters

Write the two Cartesian state coordinates as $`z=x+iy`$. The autonomous polynomial map is

```math
F(z)=\left[\lambda+i\omega+(\alpha+i\delta)|z|^2+\beta\mathrm{Re}(z^n)\right]z
     +\gamma\bar z^{n-1}+\eta\bar z^{2n-1}.
```

| Case | $`n`$ | $`\lambda`$ | $`\alpha`$ | $`\beta`$ | $`\gamma`$ | $`\omega`$ | $`\delta`$ | $`\mathrm{Re}\eta`$ | $`\mathrm{Im}\eta`$ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 13745 | 17 | 2.312147 | -2.979882 | 0.00705639 | 1.21458 | 0.1358503 | -0.3287479 | -0.2317855 | -0.1112812 |
| 205188 | 18 | 2.467358 | -3.766703 | -0.0471329 | 0.1614513 | -0.05861101 | 0.08190057 | -0.08924181 | 0.2548252 |
| 5529 | 18 | 2.286399 | -3.235571 | 0.0155766 | 0.8797143 | -0.2729181 | 0.201281 | 0.212484 | -0.324791 |
| 202080 | 20 | 2.396712 | -4.318038 | -0.04008208 | 0.4264381 | -0.02775313 | 0.01596069 | -0.1344972 | 0.2341164 |
| 1723 | 4 | 1.015272 | -0.4927374 | 0.4912107 | -0.9452272 | 0.001396406 | 0.1811445 | -0.06119981 | 0.09918609 |
| 14725 | 8 | -1.66796 | 1.365686 | -0.1338975 | 1.234103 | -0.1803585 | -0.06525403 | -0.1957832 | -0.3894122 |
| 11450 | 3 | -1.475257 | 2.268312 | 0.4538186 | 1.533973 | -0.07809093 | -0.3077413 | 0.2462413 | 0.6479135 |

The table is rounded. The [numerical receipt](numerics.json) gives each coefficient in binary64 hexadecimal, the screen seed, the seven initial states and the byte hashes. All displayed trajectories start at $`(0.247,0.193)`$, with each input decimal converted to its nearest binary64 value. The [annulus certificate](annuli.json) interprets each resulting coefficient exactly as a binary rational.

For $`\zeta^n=1`$, direct substitution gives $`F(\zeta z)=\zeta F(z)`$. Conversely, the nonzero monomial with coefficient $`\gamma`$ forces every commuting origin-centered rotation to satisfy $`\zeta^n=1`$. Complex multiplication by $`\lambda+i\omega`$ is the derivative at zero. Here $`\omega\ne0`$, so no origin-centered orthogonal reflection can commute with that derivative. This establishes the exact cyclic group for the map's origin-centered orthogonal symmetries. Its invariant components inherit the rotations. Symmetry alone does not imply that a particular orbit samples them uniformly.

Field and Golubitsky previously developed polynomial symmetric icons with the $`\lambda,\alpha,\beta,\gamma,\omega`$ terms and studied finite orbits with many rotation sectors. This search uses radial angular drift $`i\delta|z|^2z`$ and a complex second angular harmonic $`\eta\bar z^{2n-1}`$. The evidence here concerns these coefficients, explicit trapping intervals, radial permutations and recorded orbits. Historical uniqueness of the forms and families is open. [Field and Golubitsky, *Symmetry in Chaos*, Appendix A](https://epubs.siam.org/doi/abs/10.1137/1.9780898717709.appa). [Authors' paper on symmetric chaos](https://www.asc.ohio-state.edu/golubitsky.4/reprintweb-0.5/output/papers/symmetric_chaos.pdf).

## Strict radial traps and essential components

For each certified case let $`A_i=\{z:a_i\le |z|\le b_i\}`$ with the positive decimal endpoints in the certificate. Subscripts follow time around the listed cycle. The intervals are mutually disjoint. Put $`r=|z|`$, $`\psi=n\arg z`$, and $`W=e^{-i\arg z}F(z)`$. Then

```math
W=r\big[\lambda+i\omega+(\alpha+i\delta)r^2+\beta r^n\cos\psi\big]
   +\gamma r^{n-1}e^{-i\psi}+\eta r^{2n-1}e^{-2i\psi}.
```

The [Arb replay](../../../research/wow_search_20261122/symmetry/certify.py) encloses this expression on complete closed boxes in $`[a_i,b_i]\times[0,2\pi]`$. The verifier recursively bisects rational radius and turn intervals. Rounded comparisons guide partitioning only. Acceptance uses strict 100-bit ball inequalities on every leaf. For every phase it proves

```math
a_{i+1}^2+10^{-6}<|W|^2<b_{i+1}^2-10^{-6},
\qquad \mathrm{Re}W>10^{-6}.
```

The four replays use 524, 441, 557 and 1,016 leaves in table order. Therefore $`F(A_i)\subset\mathrm{int} A_{i+1}`$ with indices taken modulo $`q`$. Write $`K_i=\bigcap_{m\ge0}F^{qm}(A_i)`$. These nested nonempty compact connected images have a connected intersection. Their union $`K`$ is the attracting invariant set of the trapping neighborhood. Compactness gives $`F(K_i)=K_{i+1}`$: the forward images of $`F(A_i)`$ under $`F^q`$ lie between two successive nested images of $`A_{i+1}`$. The disjoint annuli separate the $`q`$ nonempty connected sets, proving that $`K`$ has exactly $`q`$ components. Orbits starting in the trapping neighborhood approach $`K`$ in distance, by nested compactness. Every $`K_i`$ is itself invariant under the rotations because both $`A_i`$ and $`F`$ have this symmetry.

Take a circle in $`A_i`$ wound once about zero. The positive real part of $`W`$ gives a single continuous, periodic argument for $`F(z)/z`$ on every trapping annulus. Thus applying $`F`$ to any closed loop there preserves winding number. Each $`F^{qm}(A_i)`$ contains the image of the chosen circle, still of winding one. If a path from zero to the exterior avoided $`K_i`$, truncate it beyond the outermost annulus. Its compact image would miss one of the nested compact sets and hence miss that set's winding-one loop. The winding number of this loop at points along such a path would have to change from one to zero without a crossing, which is impossible. Every component therefore separates zero from infinity. Each meets every radial ray and contains more than one point. The argument allows complicated connected sets and does not identify them as simple curves.

Direct 800-bit complex-ball propagation also certifies that the primary exact initial point enters one of the strict traps after 2, 6, 9 and 7 steps, respectively. It verifies entry for every tested start for the first three maps, and six of seven for Case 202080. The seventh five-component start exceeded the numerical radius cutoff. This proves attraction of each verified exact orbit to its designated $`K`$, without requiring that the orbit be dense in $`K`$.

## Entropy forced between the bands

Order the annuli by increasing radius. For each neighboring pair, let $`G_j`$ be the closed intervening gap, including the outer boundary of the inner annulus and the inner boundary of the outer one. These compact gaps have positive pairwise separation. Join a point of one boundary circle of a gap to the other by any continuous path contained in that gap. The images of the endpoints lie strictly inside the two target annuli specified by the permutation. The intermediate-value theorem on the radius of the image path yields a closed subpath traversing each gap between the targets. To choose such a subpath, take the first hit of the far boundary and the last preceding hit of the near boundary. Its interior stays in the chosen gap. The same argument applies in either orientation.

Starting with a radial path across $`G_j`$, apply this rule successively and restrict its parameter interval at each step. For any finite allowed word the remaining closed subinterval is nonempty and all its points follow that word. Nested compactness gives a point for every infinite one-sided word. Let $`X_+`$ be the set of points remaining in the union of the gaps at all future times. Its itinerary map is continuous by positive separation, commutes with the forward maps and is onto the one-sided path shift. Pass to $`X=\bigcap_{k\ge0} F^k(X_+)`$ for equality $`F(X)=X`$. The two graphs below have predecessors at every vertex. Any one-sided word therefore has arbitrarily long admissible pasts. Intersecting the nested compact sets $`F^k(X_+)`$ in the closed fiber of that word preserves surjectivity from $`X`$. This constructs a **factor** onto the path shift. Multiple physical points may have the same itinerary.

For Case 205188 the order $`A\to B\to C\to A`$ gives

```math
M_3=\begin{pmatrix}0&1\\1&1\end{pmatrix},
\qquad h_{\mathrm{top}}(F|_X)\ge \log((1+\sqrt5)/2)\simeq0.48121.
```

For Case 202080 write the five annuli in radial order $`A,B,C,D,E`$. Their time cycle $`A\to C\to D\to B\to E\to A`$ gives

```math
M_5=\begin{pmatrix}0&0&1&1\\0&0&0&1\\0&1&1&0\\1&0&0&0\end{pmatrix},
\quad h_{\mathrm{top}}(F|_X)\ge\log\rho(M_5)\simeq0.41401.
```

The characteristic polynomial of the second matrix is $`t^4-t^3-t^2+t-1`$, with Perron root approximately 1.5128764. All entropy bounds are per application of $`F`$. A point on a boundary of a gap maps strictly into a trapping annulus, so the invariant gap set is disjoint from the local attracting set $`K`$. Positive entropy **inside** the attracting bands remains an open question. The analogous two- and four-component gap graphs in the receipt have spectral radius one.

## Screening and independent checks

The initial Sobol screen covered 16,384 parameter rows. Its finite-time gate required a leading rate above 0.006 in every window, occupancy of at least 85 of 48×48 Cartesian bins and at least eight visits in each of 18 angular sectors. It retained 2,689 rows. Four geometric rankings gave 181 distinct previews. Cases 1723, 13745, 14725, 11450 and 5529 were selected after viewing these previews. Rechecking all retained rows for non-power-of-two radial periods found a short-lived five-band pattern in Case 15960. A 100,000-step burn-in erased the separation there, so it was kept as an external negative control.

The follow-up generated 8,192 additional settings near the scalar radius-map windows of orders three and five. It flagged 123 settings with sampled separated bands and positive leading growth. Across both screens there were 24,576 distinct parameter settings. In the first screen the numeric case is the zero-based Sobol row. In the follow-up it is 200000 plus that row. The follow-up selected two cases with separated full-angle trapping annuli and structured orbit density.

Each PNG counts two million consecutive post-transient states from one orbit in a 1,200×1,200 square with equal axes, after 100,000 discarded iterates. Color is a monotone function of visits. The complete arrays contain the raw $`x,y`$ coordinates. The same stored orbit supplies all angular sectors. For the two-, three-, four- and five-component examples, all 20 leading growth windows were positive. The strongest two-direction result was Case 13745, with second rate in the interval [0.00329, 0.00745] across its 20 windows and positive second mean at each of six additional starts. In Case 5529 each of 18 sectors received 5.55190–5.56115% of the orbit. In Case 13745 each of 17 sectors received 5.87235–5.88910%. The receipt includes the complete vectors.

For all seven candidates, [the validation program](../../../research/wow_search_20261122/symmetry/validate.py) measures growth via the analytic real-linear tangent map, checks both Jacobian columns against centered differences and replays using a separately written polar formula. It retains full-window ranges, all alternate starts, total-variation distances on a fixed 64×64 grid, a lag-1–512 late periodicity diagnostic and external trajectory hashes. Agreement of finite-time growth and visit distributions supports the numerical interpretations.

## Reproduction

From the repository root, use Python 3.11, NumPy, SciPy, Numba, Pillow and python-flint 0.9.0. The certificate replay alone needs only python-flint and the tracked JSON.

```powershell
$src = "research/wow_search_20261122/symmetry"
python "$src/certify.py" --verify
python "$src/icons.py" scout --power 14 --gallery 180
python "$src/nonpower.py"
python "$src/odd.py"
foreach ($id in 1723,13745,14725,11450,5529,205188,202080) {
    python "$src/icons.py" render --case $id --count 2000000
}
python "$src/validate.py"
```

`WOW_SYMMETRY_CACHE` overrides the default external directory `D:/Codex/chaos-attractor-cache/wow-search-20261122/symmetry`. Each `cyclic-ID.npz` contains `points` (2,000,000×2 float64), `rates` (20×2) and the complete parameter row `p`. The numerical receipt records SHA-256 for the NPZ, raw little-endian state array and PNG. [Scout and renderer](../../../research/wow_search_20261122/symmetry/icons.py). [Follow-up sampling](../../../research/wow_search_20261122/symmetry/odd.py). [Sampling-based annulus proposal](../../../research/wow_search_20261122/symmetry/annuli.py).

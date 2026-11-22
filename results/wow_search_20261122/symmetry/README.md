# Cyclic polynomial maps with nested trapping annuli

The plotted seed enters the Case 13745 attractor on step two and then alternates between two concentric bands. Two million measured steps trace 17 lobes on the inner band and a scalloped outer rim. Changing the coefficients produces three, four or five bands. Interval arithmetic proves the number of connected components and the order in which the map visits them. [Entry proof for the plotted seed](../../attractor_followup_2/annuli/README.md)

[Map and coefficients](#the-map-and-parameters) · [Proof of the bands](#strict-radial-traps-and-essential-components) · [Dynamics in the gaps](#entropy-forced-between-the-bands) · [Search](#screening-and-independent-checks) · [Reproduce](#reproduction)

The letters below label bands from the center out. Each portrait records a single orbit. [Watch 108 consecutive steps of Case 13745](alternating-orbit.webp).

| Case | Rotation order | Certified radial cycle, innermost = A | Mean numerical rates per iterate | Portrait |
| --- | ---: | --- | --- | --- |
| 13745 | 17 | A → B → A | +0.125195, +0.005589 | [Two nested crowns](cyclic-13745.png) |
| 205188 | 18 | A → B → C → A | +0.086178, +0.000386 | [Three nested rings](cyclic-205188.png) |
| 5529 | 18 | A → C → B → D → A | +0.014771, −0.001591 | [Four interleaved rings](cyclic-5529.png) |
| 202080 | 20 | A → C → D → B → E → A | +0.027593, −0.000218 | [Five nested rings](cyclic-202080.png) |

All seven portraits were drawn on **2026-11-22**. Three more settings were selected by numerical screening: [a four-cornered web (Case 1723)](cyclic-1723.png), [a hollow octagon (Case 14725)](cyclic-14725.png) and [three bent chains of loops (Case 11450)](cyclic-11450.png).

## The map and parameters

In Cartesian coordinates $`z=x+iy`$, iterate the polynomial

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

Every portrait starts at $`(0.247,0.193)`$. The table rounds the coefficients. The search and the certificate use the exact binary64 values. [Hexadecimal coefficients, seeds and starting points](numerics.json) · [Certified annuli](annuli.json)

For $`\zeta^n=1`$, substitution gives $`F(\zeta z)=\zeta F(z)`$. Conversely the nonzero $`\gamma`$ term forces $`\zeta^n=1`$ for every rotation commuting with the map. The derivative at zero multiplies by $`\lambda+i\omega`$. Since $`\omega\ne0`$, it commutes with no reflection through the origin. The exact rotation group therefore has order $`n`$.

Field and Golubitsky studied the polynomial with the $`\lambda,\alpha,\beta,\gamma,\omega`$ terms. The map here adds radial angular drift $`i\delta|z|^2z`$ and a second angular harmonic $`\eta\bar z^{2n-1}`$. [*Symmetry in Chaos*, Appendix A](https://epubs.siam.org/doi/abs/10.1137/1.9780898717709.appa) · [Field and Golubitsky on symmetric chaos](https://www.asc.ohio-state.edu/golubitsky.4/reprintweb-0.5/output/papers/symmetric_chaos.pdf)

## Strict radial traps and essential components

For each certified case let $`A_i=\{z:a_i\le |z|\le b_i\}`$ with the positive decimal endpoints in the certificate. Subscripts follow time around the listed cycle. The intervals are mutually disjoint. Put $`r=|z|`$, $`\psi=n\arg z`$, and $`W=e^{-i\arg z}F(z)`$. Then

```math
W=r\big[\lambda+i\omega+(\alpha+i\delta)r^2+\beta r^n\cos\psi\big]
   +\gamma r^{n-1}e^{-i\psi}+\eta r^{2n-1}e^{-2i\psi}.
```

The verifier covers $`[a_i,b_i]\times[0,2\pi]`$ with closed boxes. Rounded estimates guide subdivision. On each final box, 100-bit Arb interval arithmetic proves both strict inequalities. [Verifier](../../../research/wow_search_20261122/symmetry/certify.py)

```math
a_{i+1}^2+10^{-6}\lt |W|^2\lt b_{i+1}^2-10^{-6},
\qquad \mathrm{Re}W\gt 10^{-6}.
```

**Component count.** The four cases required 524, 441, 557 and 1,016 boxes in table order. The first inequality gives $`F(A_i)\subset\mathrm{int} A_{i+1}`$, cyclically modulo $`q`$. Each $`K_i=\bigcap_{m\ge0}F^{qm}(A_i)`$ is a nested intersection of nonempty compact connected sets. The successive images of $`F(A_i)`$ under $`F^q`$ lie between successive images of $`A_{i+1}`$, so $`F(K_i)=K_{i+1}`$. Their union $`K`$ attracts the trapping neighborhood by nested compactness. The disjoint annuli separate its $`q`$ connected components. Each inherits the rotational symmetry.

**Each component surrounds zero.** The second inequality gives a continuous periodic argument for $`F(z)/z`$ on every trapping annulus. Applying $`F`$ therefore preserves the winding number of a closed loop there. Every $`F^{qm}(A_i)`$ contains the image of a circle wound once around zero. If a path reached the exterior from zero while avoiding $`K_i`$, its compact portion up to the exterior would miss one of the nested sets. It would then avoid this set's winding-one loop even as the loop's index along the path changed from one to zero, a contradiction. Thus each component meets every radial ray. Whether it is a simple curve remains open.

**Entry from the plotted start.** At 800-bit precision the starting point enters a trapping annulus after 2, 6, 9 and 7 steps, respectively. All seven starts were certified for the first three maps and six for Case 202080. The seventh exceeded the numerical radius cutoff. Certification establishes attraction to the component cycle. Density of a single orbit in that cycle remains open.

## Entropy forced between the bands

For the three-band map, every infinite itinerary through the two gaps that avoids consecutive visits to the inner gap occurs. For the five-band map the allowed itineraries have the four-state graph below. These claims concern points that stay **between** the attracting bands.

**Crossing a gap.** Let $`G_j`$ be the closed gap between neighboring annuli, including its two boundary circles. Distinct gaps have positive separation. Join the boundaries of one gap by a path in it. Its endpoints map into the interiors of their respective target annuli. By continuity of radius the image crosses every intervening gap, in either orientation. For each crossing, take the first hit of the far boundary and the last preceding hit of the near boundary. The subpath stays in the chosen gap.

**Coding a trajectory.** Start with a radial path across $`G_j`$ and restrict its parameter interval at each crossing. Every finite admissible word leaves a nonempty closed interval. Nested compactness realizes each infinite one-sided word. Let $`X_+`$ contain the points whose forward gap labels follow the graph. This is compact because the finite union of allowed closed transition conditions is closed at each time. Positive gap separation makes the itinerary map continuous, and it intertwines $`F`$ with the one-sided shift.

To obtain an invariant domain, set $`X=\bigcap_{k\ge0}F^k(X_+)`$. Each vertex has a predecessor, so any word has arbitrarily long admissible pasts. Its closed fiber meets every nested $`F^k(X_+)`$ and hence meets $`X`$. The resulting itinerary map is onto, with $`F(X)=X`$. The entropy of the symbolic shift gives a lower bound for the entropy on $`X`$.

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

The second matrix has characteristic polynomial $`t^4-t^3-t^2+t-1`$ and Perron root 1.5128764 to the displayed precision. The bounds are per iteration. A gap-boundary point moves inside a trapping annulus, so $`X`$ is disjoint from the local attractor $`K`$. For Case 13745 a [separate inverse-branch proof](../../low_hanging_20261122/annuli/README.md) gives positive entropy inside the local attractor and proves that each component contains an open set. Positive entropy within the other three cases’ attracting bands remains open. The two- and four-band gap graphs have spectral radius one.

## Screening and independent checks

The first Sobol search tested 16,384 parameter sets. To pass, an orbit needed a leading growth rate above 0.006 in every window, at least 85 occupied bins on a 48×48 grid, and at least eight visits in each of 18 angular sectors. Of 2,689 survivors, four geometric rankings yielded 181 previews. Visual inspection selected Cases 1723, 13745, 14725, 11450 and 5529. In Case 15960 an apparent five-band pattern disappeared after 100,000 steps.

A second search sampled 8,192 settings near period-three and period-five windows of the scalar radius map. Of 123 candidates with separated bands and positive leading growth, Cases 205188 and 202080 admitted full-annulus certificates. Case numbers are zero-based Sobol row numbers, with 200000 added for this second search.

Each 1,200×1,200 portrait counts two million consecutive visits after discarding the first 100,000. Brighter pixels received more visits. The 17 sectors in Case 13745 each received 5.87235–5.88910% of the orbit. The 18 sectors in Case 5529 each received 5.55190–5.56115%.

The leading growth rate was positive in all 20 windows for each certified map. For Case 13745 the second rate ranged from 0.00329 to 0.00745, with positive means at six additional starts. For Case 205188 the second rate crossed zero across windows. These finite measurements do not settle chaos on the attracting components.

The validator compares the analytic tangent map with centered differences and replays each orbit in polar coordinates. It records the growth windows, independent starts, 64×64 occupancy comparisons and a periodicity check at lags 1–512. [Validation code](../../../research/wow_search_20261122/symmetry/validate.py) · [Results and array hashes](numerics.json)

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

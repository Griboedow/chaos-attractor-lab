# Fifteen basins on one Newton boundary

Six starting points can converge to six different roots of a single polynomial. Nine more can converge to nine phases of one cycle when sampled every ninth step. For Case 6902, all fifteen outcomes occur in a single disk of radius $`3\times10^{-32}`$. An interval certificate gives an open box of starting points for each outcome. The same fifteen basins meet at every point of the Julia set. [Interval certificate](certificate.json)

![Fifteen Newton basins in an ordinary window and in a microscopic window centered on a pole](fifteen-newton-basins.png)

<!-- drawing-date -->
*Drawn: **2026-11-22**.*
<!-- /drawing-date -->

Case 6902 identifies the coefficient isolated below. The left window covers $`[-1.85,1.85]\times[-1.20,1.20]`$ in the starting plane. The right window has native width $`2\times10^{-30}`$, centered on the marked pole. Its central ring encloses the fifteen certified sets. Cool colors identify roots. Warm colors distinguish the phases of the cycle. Darker shading within each color indicates slower capture relative to the earliest capture in that panel.

## The map and its outcomes

```math
p_a(z)=z^6-az+a,\qquad N_a(z)=z-\frac{p_a(z)}{p_a'(z)}=\frac{5z^6-a}{6z^5-a},\qquad G=N_a^9.
```

The displayed coefficient is $`a\simeq1.277961358032196+0.074550069827001i`$. Its exact definition is the unique zero of $`N_a^9(0)`$ in a complex disk of radius $`10^{-142}`$ around the dyadic center in the certificate. All nine iterates have distinct positions with pairwise separation greater than 0.003642. All six roots have pairwise separation greater than 0.6017. The coefficient is complex, giving the pictured basin pattern its asymmetry.

For $`N_a`$ there are seven attracting cycles. Each simple root is a superattracting fixed point. The remaining critical point at zero has multiplicity four and belongs to the superattracting nine-cycle. This accounts for all ten critical points with multiplicity. In particular the map is hyperbolic. For $`G`$ the same cycle supplies nine distinct attracting fixed points, giving fifteen in total. The common boundary is a forward repeller and can be sampled by inverse iteration.

Write $`B_j`$ for the **full basin** under $`G`$ of one of the fifteen fixed points. Each $`B_j`$ is open, completely invariant, and a union of Fatou components. Any neighborhood of a Julia point that missed $`B_j`$ would have all its forward images miss the whole open basin. Those images would omit three fixed values on the sphere, hence form a normal family by Montel's theorem. This contradicts the definition of a Julia point. Conversely a Fatou component is open and lies in one basin, so it contains no basin-boundary point. Thus every $`\partial B_j`$ equals $`J(G)=J(N_a)`$, the Wada property for these full basins. A pointwise statement about boundaries of individual immediate components would require a different argument. [Milnor's lecture notes on one-variable dynamics](https://arxiv.org/abs/math/9201272)

## Fifteen choices at one pole

The selected pole is $`\xi=(a/6)^{1/5}\simeq0.734161053096716+0.008556170353017i`$ on the principal branch. It maps to the repelling fixed point at infinity. A 656-bit Arb calculation proves contraction for the coefficient equation, for the six root equations, and on capture disks around the six roots and zero. On the zero-centered disk of radius $`10^{-7}`$, the modulus of the derivative of $`N_a^9`$ is below $`7.18\times10^{-24}`$.

Fifteen starting boxes of coordinate half-width $`10^{-101}`$ lie inside $`|z-\xi|<3\times10^{-32}`$. Direct interval propagation places each entire box in its assigned capture disk after 386–394 ordinary Newton steps. A second evaluator using 160-digit `Decimal` arithmetic and the original polynomial/derivative form agrees on every center and capture time. The boxes were deliberately constructed by inverse branches. Their counts do not estimate area fractions. [Witness boxes and exact dyadics](certificate.json) · [Independent replay](independent-decimal.json)

### A zoom that permutes nine colors

Near infinity, in the coordinate $`1/z`$, the Newton map has multiplier $`6/5`$. Let $`L`$ be its local inverse branch fixing infinity. At $`\xi`$ the Newton map has a simple pole, so its restriction $`H`$ there is locally invertible on the sphere. The local change of coordinates

```math
T=H^{-1}\circ L\circ H,\qquad T(\xi)=\xi,\quad T'(\xi)=5/6,\qquad N_a^2\circ T=N_a
```

preserves the six root labels and moves a point in phase $`j`$ into phase $`j-1\pmod9`$. Applying it nine times restores all colors. In the original plane, simple Euclidean contraction by $`5/6`$ is its first-order approximation.

We compared 204 identical offsets between a window centered at the pole and a copy shrunk by $`5/6`$ using direct Arb iteration. The expected color permutation held for 136, 186, 204, and 204 pairs at respective initial half-widths $`10^{-1}, 10^{-3}, 10^{-8}, 10^{-30}`$. All 204 shrunken starts at the finest scale reached their capture set one step later. The sample combined 144 uniform pixel centers with four per outcome in the finest window. [Paired interval results](scale-control.json)

## Independent geometric checks

The full-resolution left image has 1,420,800 classified starts. The right image has 800,000. For the microscopic raster, the first step uses the principal part of the Laurent expansion near the pole. Iterates above modulus $`10^5`$ are advanced using the leading term $`(5/6)z`$ before switching to the rational formula. On 299 selected pixel centers, direct interval iteration agreed on both color and capture time. These include a uniform random sample and representatives of every outcome. [Forward screen](inverse-check.json) · [Microscopic screen](microscopic-check.json)

Independently, inverse iteration solves $`5z^6-6wz^5+a(w-1)=0`$ for all six branches over $`w`$. Starting from the five finite poles gives 1,399,680 numerical preimages that first reach infinity at depth eight. Deeper exploration retaining one point per occupied voxel of side 0.004 on the unit sphere produces 99,423 representatives through depth 23. The Julia set is the closure of the backward orbit of the nonexceptional point at infinity. Both saved clouds cover the sphere without a planar crop. For $`z=x+iy`$, the sphere coordinates are $`(2x,2y,x^2+y^2-1)/(1+x^2+y^2)`$.

All 67,760 retained inverse representatives inside the one-pixel interior margin of the left window lie within two pixels of a forward color interface. Conversely, 83.68% of 160,000 uniformly sampled interface pixels lie within two pixels of a retained inverse point. Median distance is 0.834 pixels, and the 99th percentile is 22.62 pixels. This asymmetry records the finite inverse coverage. The unstratified depth-eight cloud alone covers 53.22% of the same interface sample within two pixels. [Inverse comparison](inverse-check.json) · [Depth-23 control](stratified-inverse.json)

## Search and historical context

A scrambled Sobol search seeded 8,192 coefficient solves in each of 21 degree/period cells, with degrees 4–6 and requested periods 3–9. Of 172,032 starts, 114,778 met the numerical solve criterion. The screen grouped conjugate coefficients by six-decimal rounding, inspected up to 50 clusters per cell, and tested 801 clusters in total. Fifty-four passed the coarse spatial and visibility filters. Lower periods and repeated parameters can occur across cells. We selected this sextic for its nine visible phase regions and certified its exact period separately. [Search counts, array shapes and hashes](manifest.json)

Campbell and Collins previously proved that this polynomial family supports any prescribed period for every degree at least three. Their parameterization is obtained by setting $`a=(d-1)c^d`$. The contribution here is the particular complex coefficient and the linked numerical and interval experiments on one basin boundary. The period, family, Wada mechanism and inverse iteration have established precedents. [Campbell and Collins, *Specifying attracting cycles for Newton maps of polynomials*](https://arxiv.org/abs/1208.5687) · [Milnor, *Dynamics in One Complex Variable*](https://arxiv.org/abs/math/9201272)

## Reproduction

Run the following PowerShell commands from the repository root. `WOW_BOUNDARY_CACHE` names the directory for raw arrays. The recorded environment used Python 3.12.14, NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.0, Pillow 13.0.0 and python-flint 0.9.0. The certificate replay works directly from the checked-in JSON. Search results and the raw meshes can be regenerated independently of that replay.

```powershell
$src = 'research/wow_search_20261122/boundary'
$out = 'results/wow_search_20261122/boundary'
$env:WOW_BOUNDARY_CACHE = 'D:/Codex/chaos-attractor-cache/wow-search-20261122/boundary'
$env:NUMBA_CACHE_DIR = "$env:WOW_BOUNDARY_CACHE/numba"
$env:OPENBLAS_NUM_THREADS = '1'
New-Item -ItemType Directory -Force $env:WOW_BOUNDARY_CACHE | Out-Null
python "$src/certificate.py" --verify "$out/certificate.json"
python "$src/search.py" --levels 13 --max-seen 50
python "$src/inverse.py" --depth 8
python "$src/stratified.py"
python "$src/certificate.py"
python "$src/decimal_control.py"
python "$src/microscope.py" --width 1000 --height 800
python "$src/scale_control.py"
python "$src/figure.py"
python "$src/receipt.py"
```

The primary catalog cloud is `case-6902-stratified-inverse-sphere.npz`, with `points[99423,3]`, `plane[99423,2]` and `first_infinity_depth[99423]`. Separate archives retain all depth-eight points in native and sphere coordinates. Standard and microscopic label grids retain the capture iterations. Raw paths are relative to the parent of `WOW_BOUNDARY_CACHE` in the manifest. Floating-point branch ordering and image fonts can vary by platform. [Source files](../../../research/wow_search_20261122/boundary) · [Machine-readable receipt](manifest.json)

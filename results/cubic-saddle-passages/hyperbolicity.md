# A horseshoe on three basin boundaries

In the [forced cubic at gain 0.34](../shape_search_20261121/frontiers/cubic/#a-single-gradient-field), the [six-phase diagram](README.md) codes a compact, uniformly hyperbolic subset of the common boundary. **Every two-sided route picks out exactly one point in this subset at phase zero.** It contains exactly 18 fixed saddles and 152 distinct saddle cycles of least period two. Its entropy is $`\log(9+4\sqrt5)`$ per forcing period. The full set of orbits visiting the larger rectangles may contain additional points.

<details>
<summary>How to identify one point per route</summary>

In the diagram's local coordinates, prescribe the longitudinal coordinate $`X`$ at the **end** of a passage and the transverse coordinate $`Y`$ at its **start**. For each arrow, the checker constructs a smooth map $`C(X,Y)=(x_{\rm start},y_{\rm end})`$ into the interior of the square $`[-0.14,0.14]^2`$. Every selected map satisfies $`\lVert DC\rVert_\infty<0.635195`$ throughout that square.

For a fixed route, update each $`x_i`$ from $`C_i(x_{i+1},y_i)`$ and each $`y_i`$ from $`C_{i-1}(x_i,y_{i-1})`$. This simultaneous update contracts in the supremum norm. Its unique fixed sequence follows the ODE. After $`n`$ updates each coordinate depends only on a window of radius $`n`$, so the solution varies continuously with the itinerary. Different routes visit disjoint rectangles at some phase. This gives a conjugacy. The same mixed derivative bound maps the horizontal tangent cone strictly into itself and expands it by at least $`1/0.635195`$ in longitudinal coordinate. Applying it backwards gives the vertical cone estimate. [Gonchenko and Li on the Shilnikov cross-map criterion, §3.2](https://doi.org/10.1134/S1560354710020061)

The constant right-site code selects $`+1`$, whose [stable manifold lies on all three boundaries](../cubic-saddle-basins/). Routes eventually following this code are dense and converge forward to $`+1`$. Closedness puts the entire horseshoe on the common boundary.

</details>

<details>
<summary>What the interval computation checks</summary>

Each passage has 32 shooting segments. On each of 256 overlapping parameter patches, an interval Newton contraction encloses all 33 nodes for every pair of prescribed endpoints. Taylor defects and a second-variation bound enclose the shooting Jacobian throughout each node box. Implicit differentiation bounds both row sums of the mixed map. On each adjacent pair of patches the enclosed center solution of one fits in the uniqueness box of the other. Local uniqueness and connectedness then identify their branches on the overlap. Nine arrows require computation and half-period reflection supplies the other nine.

The checker also transports the [18 archived periodic saddle enclosures](../cubic-forced-orbits/) through each passage and checks every shooting node against the selected branch's uniqueness box.

[Interval bounds by arrow](hyperbolicity.json) · [Checker](../../research/shape_search_20261121/frontiers/cubic/boundary_orbits_graph.py)

</details>

```sh
python research/shape_search_20261121/frontiers/cubic/boundary_orbits_graph.py --workers 7
```

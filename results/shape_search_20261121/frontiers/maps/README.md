# Positive entropy on ten connected components

Case 2300491 has a localized attracting set A with exactly ten nondegenerate connected components. For every fixed coefficient vector in the closed ten-dimensional box of radius $`2^{-23}`$ around its recorded coefficients, a compact invariant subset of the twentieth iterate maps continuously onto the full two-sided shift on two symbols. The uniform conclusions include

```math
h_{\mathrm{top}}(F|_A)\geq\frac{\log 2}{20}=0.034657359\ldots,
\qquad 1.018256<\dim_H A\leq1.6.
```

Every component has topological covering dimension one and the same Hausdorff dimension as A. The [entropy certificate](entropy-certificate.json) gives the native affine frame, exact model coefficients, intermediate boxes, complete dyadic partitions, strict margins and face controls. The [dimension certificate](dimension-certificate.json) gives the native phase hulls and derivative enclosures used for the dimension estimates. Decimal model and chart coefficients are interpreted as the binary64 values obtained by parsing the JSON. These results hold for each fixed parameter choice. Uniqueness of the symbolic code's realization, hyperbolicity and a positive asymptotic exponent for the canonical orbit remain open.

## What is certified

Let G be the tenth iterate of the three-axis twist map, T the closed trapping box, and H the twentieth iterate. In its normalized frame T is the cube $`[-1,1]^3`$. The directed intervals and 144-bit Arb calculations prove $`G(T)\subset\operatorname{int}T`$ and $`F^j(T)\cap T=\varnothing`$ for all nine earlier positive iterates. Thus the nested intersection $`A_0=\bigcap_{k\ge0}G^k(T)`$ is connected and attracts a neighborhood under G. The map is invertible because each twist can be reversed at its preserved axis coordinate and the affine ratio is positive. Its first ten phase images of A₀ are mutually disjoint. They form the ten components of A.

Write u for the third normalized coordinate. The two closed slabs occupy longitudinal intervals $`[-0.5625,-0.25]`$ and $`[-0.1875,0.28125]`$, each across the full transverse square. For each of their four faces, the certificate encloses its first G-image in the corresponding intermediate box below. It places the G-image of that entire intermediate box strictly beyond the indicated outer longitudinal face.

| Face u | Intermediate coordinate 1 | Intermediate coordinate 2 | Intermediate u | Final u |
|---:|---|---|---|---|
| −0.5625 | [−0.200, −0.147] | [−0.355, −0.283] | [−0.417, −0.308] | < −0.5625 |
| −0.2500 | [−0.018, 0.082] | [−0.096, 0.094] | [−0.940, −0.787] | > 0.28125 |
| −0.1875 | [−0.010, 0.101] | [−0.107, 0.097] | [−0.961, −0.797] | > 0.28125 |
| 0.28125 | [−0.122, −0.069] | [−0.723, −0.590] | [−0.461, −0.296] | < −0.5625 |

The two entire three-dimensional source slabs lie in T. Since G maps all of T strictly into its interior, H maps both complete slabs there. In particular every returning point satisfies both strict transverse inequalities. The largest tree for these conclusions has 901 nodes and covers the complete trapping box, including its boundary.

Any continuous path joining the two longitudinal faces of one slab maps to a path crossing either target slab. Orient the image to start below the target. Take its first visit to the upper target coordinate and the last preceding visit to the lower one. The subpath between them lies in the target slab by the transverse bounds. Restriction and iteration realize every finite prescribed itinerary. To realize a two-sided infinite itinerary, realize its symbols from time −n through n and retain the point at time zero. The compatible closed sets are nested and nonempty in a compact slab. Their intersection realizes the entire itinerary. Invertibility and trapping put the resulting compact H-invariant set in A₀. The slab labels define a continuous surjection onto the binary shift because the slabs are compact and separated. The factor inequality and entropy of an iterate give the stated bound for F. Multiple points may share an itinerary.

## Interval calculations

The complete parameter block is carried through the proof. Each cell is propagated around an exactly represented nominal center. At step i, enclose the native derivative Jᵢ over the convex hull of the center and cell, choose a point matrix Cᵢ, and include the single-step point residual bᵢ from the full coefficient block. For the deviation e from the center the mean-value theorem gives an enclosure

```math
e_i\in C_i e_{i-1}+(J_i-C_i)e_{i-1}+b_i .
```

The verifier transports each remainder separately through subsequent center matrices and intersects with the natural interval enclosure. Arb regenerates all nine complete partitions. A separate directed binary64 implementation uses rational sine and cosine bounds from Machin's identity and Taylor expansions, then reads each tree to its end. The four face-center controls also disprove the opposite crossing orientations.

## Hausdorff dimension

For each fixed longitudinal coordinate in the first slab, and for every infinite past, the crossing-path argument gives a point in that exact slice with that past. At finite depth end with a path spanning the first slab, and then intersect that path with the prescribed plane. Compactness preserves the exact slice and all past symbols simultaneously. Both slabs have a certified positive native gap greater than 0.001920643.

For each phase a complete dyadic covering of T encloses its forward image in a convex native box. A line segment joining any two points of that phase stays in this box. Their inverse difference is their initial difference multiplied by the average inverse Jacobian along the segment. Each of the ten average matrices is therefore contained entrywise in the derivative hull over its entire convex phase box. The two derivative implementations each cover 4,096 subboxes per phase. Compose the interval matrices in inverse order and bound the result in a fixed weighted maximum norm at phase zero. Squaring this bound gives a Lipschitz constant L for the inverse of H between any two points of A₀. The common Arb phase hulls serve both inverse-derivative implementations.

Different pasts first disagreeing k positions back give points separated by at least a constant times $`L^{-k}`$ in native coordinates. Equivalence of the fixed weighted norm and native distance contributes a single constant independent of k. Thus the codes represented in any sufficiently small set of diameter r have binary Bernoulli outer mass bounded by $`C r^{\log(2)/\log L}`$, with the same C for all of the exact slices. Every slice realizes all codes, so every slice has a uniform positive Hausdorff content at that exponent. The longitudinal projection of any native covering ball has length at most a constant times its diameter. Integrating the slice content inequality over the positive-length first interval gives the dimension lower bound $`1+\log(2)/\log L`$. No selection of one realizing point per past is required.

For the upper bound, every source point of T satisfies a uniform strict singular-value inequality for DG. For $`d=1+q`$, the code bounds $`\sigma_1\sigma_2^q`$ by $`\lVert DG\rVert_F^{1-q}\lVert\operatorname{cof}DG\rVert_F^q`$ on each dyadic cell. The two implementations propagate native tangent and exterior-square matrices through all ten iterates, with outward rounding and bounds from the direct minors. The uniform singular-value covering bound for the compact invariant A₀ yields $`\dim_H A_0\le d`$. The finite phase maps are bi-Lipschitz on compact neighborhoods, so all components have that dimension bound. Connectedness and the nontrivial symbolic subsystem ensure topological dimension at least one; the Hausdorff upper bound below two gives the opposite bound for topological dimension.

## Reproduction

From the repository root, use Python 3.11 with NumPy, SciPy, Numba and python-flint. The single command below verifies the tracked proof, without the external point-cloud cache. It regenerates the Arb partitions, replays the recorded partitions in directed intervals, repeats the upper-dimension verification, checks the negative controls and verifies the evidence hashes. The [manifest](manifest.json) records the selected bounds.

```powershell
python research/shape_search_20261121/frontiers/maps/packet.py verify results/shape_search_20261121/frontiers/maps
```

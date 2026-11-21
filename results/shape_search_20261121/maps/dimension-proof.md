# A noninteger Hausdorff dimension for the five-component attracting set

Two independent interval implementations establish the conservative bounds $`1.039<\dim_H A\le1.45`$ throughout the closed parameter block of radius $`2^{-24}`$ around each of the ten coefficients. Each component has topological covering dimension one. The canonical start enters the trapping box at step six for every parameter tuple in this block.

Let $`F`$ be Case 2503155 and let $`G=F^5`$. The closed return box $`T`$ has $`G(T)`$ strictly inside its interior, and $`T`$ misses its first four forward images. Write $`A_0=\bigcap_{n\ge0}G^n(T)`$ and $`A=\bigcup_{i=0}^4F^i(A_0)`$. The constant parameters may vary independently throughout the block specified by the two uniform certificates below. The affine frame and all boxes stay fixed throughout that block. The conclusion concerns the entire localized attracting set.

The trapping test was checked using [directed intervals](certificates/2503155-cycle-parameter.json) and [Arb](certificates/2503155-cycle-parameter-arb.json). The symbolic test was checked using [directed intervals](certificates/2503155-symbolic-parameter.json) and [Arb](certificates/2503155-symbolic-parameter-arb.json). The [exact affine comparison](certificates/2503155-containment.json) puts the symbolic frame inside the return frame. The crossing and avoidance claims use strict inequalities on closed input boxes.

## Every slice realizes every past

Use native Euclidean distance in three-space. In the symbolic frame write $`u`$ for the third normalized coordinate, an affine linear functional with constant gradient $`\ell`$. Let $`S_0`$ and $`S_1`$ be the full transverse slabs defined by the [four faces](certificates/2503155-faces.json). Their longitudinal intervals are $`[f_0,f_1]`$ and $`[f_2,f_3]`$, with $`f_0<f_1<f_2<f_3`$. Their transverse coordinates range over $`[-1,1]^2`$. Consequently their native separation $`\delta`$ satisfies $`\delta\ge(f_2-f_1)/\|\ell\|>0`$.

For either source slab, $`G`$ takes its two longitudinal faces strictly past opposite outer faces $`f_0`$ and $`f_3`$. For every point of the entire source slab whose output has $`u\in[f_0,f_3]`$, both output transverse coordinates lie strictly in $`(-1,1)`$. Take any path from one longitudinal face of a source slab to the other. Orient its image to start below $`f_0`$. For a target interval $`[a,b]`$, take the first time its longitudinal coordinate reaches $`b`$, then the last time before that when it equals $`a`$. The image segment between those times stays in $`[a,b]`$. The certified transverse bounds put this entire segment in the target slab. The image segment crosses that slab from face to face. Each restriction retains the earlier prescribed slabs in the path's past.

Fix any exact $`t\in[f_0,f_1]`$ and any infinite binary past $`(b_1,b_2,\ldots)`$, with $`b_1`$ the most recent previous symbol. For a finite depth $`N`$, start with a longitudinal crossing in $`S_{b_N}`$. Apply the preceding construction in the order $`S_{b_{N-1}},\ldots,S_{b_1},S_0`$. The final crossing meets the plane $`u=t`$ by the intermediate value theorem. Invertibility gives this intersection point its specified finite past. In particular, each of the compact sets

```math
K_N(t,b)=S_0\cap\{u=t\}\cap\bigcap_{j=1}^N G^j(S_{b_j})
```

is nonempty. They decrease with $`N`$, so their intersection contains a point $`x`$. Its full backward return orbit remains in the supplied slabs inside $`T`$. For every $`n`$, the equality $`x=G^n(G^{-n}x)`$ puts $`x`$ in $`G^n(T)`$, hence in $`A_0`$. This constructs a point in the exact chosen slice for every independently chosen infinite past, including the extreme slices $`t=f_0,f_1`$. Multiple points may realize one past.

## Native inverse Lipschitz estimate

The first slab extends from $`-0.44791666666666674`$ to $`0.03125`$, and the second from $`0.05208333333333326`$ to $`0.6354166666666665`$. The gradient norm is below 10.25, giving native separation greater than 0.00203.

Every point of $`A`$ has norm at most one. Indeed, $`F`$ maps $`A`$ onto itself. The norm after one step is bounded by $`r\|x\|+1-r`$, since the three rotations preserve it. Taking the maximum over compact $`A`$ gives the claim for $`0<r<1`$. The native unit ball is convex.

Invert a step by reversing its three spatial rotations and then undoing its affine contraction. Each reverse rotation preserves the native unit ball. In a participating plane its derivative is an orthogonal rotation composed with a shear of size at most $`|k+3hw^2|\sqrt{1-w^2}`$ at unit-ball points, where $`w`$ is the unchanged perpendicular coordinate. To see this, differentiate the rotated two-vector with respect to $`w`$. The resulting column has length $`|k+3hw^2|`$ times the length of that two-vector. The largest singular value of the resulting three-dimensional shear with column length $`v`$ is $`(\sqrt{4+v^2}+v)/2`$. The inverse affine step has norm $`1/r`$.

The [lower-bound calculation](certificates/2503155-dimension-lower.json) bounds $`|k|`$ and $`3|h|`$ for each of the three twist stages over the whole parameter block. For their nonnegative upper bounds $`K,H`$, the squared shear is bounded by the maximum of $`(K+Hy)^2(1-y)`$ on $`[0,1]`$. It suffices to test the endpoints and the internal critical point $`(2H-K)/(3H)`$ when present. These comparisons use exact rationals. Arb encloses the square roots, products and logarithms. Multiplying the three stage bounds by the reciprocal lower bound for $`r`$ gives a one-step inverse Lipschitz bound $`L_1`$ on the native unit ball. The line segment joining any two points of $`A`$ lies in this convex domain. Both of their actual inverse images stay in $`A`$, so repeating the one-step estimate gives

```math
\|G^{-j}x-G^{-j}y\|\le L^j\|x-y\|,\qquad x,y\in A_0,\quad L=L_1^5.
```

## Uniform positive content and integration over slices

The global inverse estimate alone gives $`\dim_H A>1.022934`$. A sharper estimate encloses each compact phase $`A_i`$ in a native axis-aligned box $`H_i`$. Starting with the full closed $`T`$, propagate every cell of a dyadic covering to step $`i`$ and take the coordinate hull of the results. This covers $`F^i(T)`$ and hence $`A_i`$. Given two points in $`A_i`$, the entire straight segment joining them lies in the convex $`H_i`$. Their inverse difference equals their initial difference multiplied by the integral average of $`DF^{-1}`$ along this segment. Every entry of that average belongs to the corresponding closed interval of an enclosure of $`DF^{-1}(H_i)`$. A separate dyadic partition of $`H_i`$ bounds all nine entries.

Starting in $`A_0`$, five inverse steps depart from phases $`0,4,3,2,1`$ in that order. Multiply the five interval matrices in application order and bound the resulting matrix's operator norm by its Frobenius norm. The result is a Lipschitz bound for $`G^{-1}`$ between any two points of $`A_0`$. Composing it $`j`$ times gives the same inverse-distance inequality above. All bounds include the full parameter block. [Arb phase enclosures](certificates/2503155-dimension-phase-arb.json), [independent directed-interval phase enclosures](certificates/2503155-dimension-phase-directed.json).

Each implementation used 262,144 dyadic subboxes per region. The native inverse-return bounds are below 34,620,000 in Arb and below 29,511,000 in directed binary64 intervals. The corresponding strict rational Hausdorff lower bounds are 1.039928 and 1.040298, respectively. Both computations establish the conservative lower bound stated at the start of this proof.

Let $`E_t`$ consist of points in $`S_0\cap\{u=t\}`$ with a full past through the two slabs. It is a compact subset of $`A_0`$. Invertibility and the positive gap between slabs give each such point a unique past. The preceding construction makes the coding onto all one-sided binary sequences. If two pasts first differ at index $`j`$, their inverse images at that index lie at least $`\delta`$ apart. Thus their two current native points lie at least $`\delta L^{-j}`$ apart.

Put $`s=\log 2/\log L`$. Under fair Bernoulli measure, the codes from any native ball of radius $`a`$ have mass at most $`2(2a/\delta)^s`$. Indeed, if $`2a<\delta L^{-n}`$, all codes in that ball share their first $`n`$ symbols, a cylinder of mass $`2^{-n}`$. Choose the largest such $`n`$. Then $`2a\ge\delta L^{-(n+1)}`$, giving the claimed bound. For a large ball the same bound is trivial. Cylinder containment also gives the bound for outer measure without any measurability assumption on the codes of a ball.

For any countable ball cover of $`A_0`$ with native radii $`a_i`$, the subcollection of balls meeting the plane $`u=t`$ covers $`E_t`$. Their codes cover all pasts and have total Bernoulli outer measure at least one. Therefore

```math
\sum_{i:\,B_i\cap\{u=t\}\ne\varnothing} a_i^s\ge c,
\qquad c=\tfrac12(\delta/2)^s>0
```

for every $`t\in[f_0,f_1]`$, with the same constant. A ball of radius $`a_i`$ meets planes over an interval of $`t`$ of length at most $`2\|\ell\|a_i`$. Integrating the preceding inequality over the target interval and summing nonnegative terms gives

```math
\sum_i a_i^{1+s}\ge \frac{c(f_1-f_0)}{2\|\ell\|}>0.
```

This holds on arbitrarily fine covers. The $`(1+s)`$-dimensional Hausdorff measure of $`A_0`$ is positive, and $`\dim_H A_0\ge1+s`$. Each $`F^i`$ is a smooth diffeomorphism and is locally bi-Lipschitz on compact sets. Its five compact images have the same Hausdorff dimension, which is also the dimension of their finite union $`A`$.

## Uniform upper bound

For $`1<d<2`$ write $`q=d-1`$. At each point $`x`$ of the entire closed return box, let $`J=DG(x)`$ and $`C=\operatorname{cof}J`$. If the singular values of $`J`$ are $`\sigma_1\ge\sigma_2\ge\sigma_3`$, the largest singular value of $`C`$ is $`\sigma_1\sigma_2`$. The [Arb certificate](certificates/2503155-dimension-arb.json) verifies a uniform strict upper bound below one for

```math
\sigma_1\sigma_2^q
\le \|J\|_F^{1-q}\|C\|_F^q.
```

The subscript $`F`$ denotes the Frobenius norm. The checker propagates a native tangent and its exterior square through the return, intersects the latter with direct two-by-two minors, and compares the 40th power of the displayed bound with one in every leaf of a complete dyadic partition. The exponent is $`d=29/20`$, so $`q=9/20`$. The Arb computation needed 21,567 nodes and gave a largest accepted upper bound below 0.99885. A [separate directed-interval computation](certificates/2503155-dimension-directed.json) bounds cofactors using products of single-step minors and accepted 44,586 leaves with a largest bound below 0.99984. Both upper-bound tests allow a parameter block of radius $`2^{-22}`$, which contains the stated block. The Douady–Oesterlé singular-value criterion gives $`\dim_H A_0\le d`$ for the compact invariant $`G(A_0)=A_0`$. [Original theorem](https://mathscinet.ams.org/mathscinet-getitem?mr=585918), [generalized singular-value criterion](https://dml.cz/manakin/handle/10338.dmlcz/703054).

The strict lower bound above one and the upper bound below two place the Hausdorff dimension strictly between consecutive integers. Each connected component has at least two points because the disjoint symbolic slabs each contain a periodic point. Its topological covering dimension is therefore at least one and is at most its Hausdorff dimension. Its topological covering dimension is exactly one. The same conclusions hold for the union of its five components.

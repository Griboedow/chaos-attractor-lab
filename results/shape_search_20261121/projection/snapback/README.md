# Certified positive-entropy subsystems

Four parameter sets have a compact invariant subsystem of positive topological entropy. Each packet records two continuous inverse branches with disjoint images. The bounds are per original map update.

| Forward portrait | Surface | Least repelling period | Common return | Certified lower bound |
|---|---|---:|---:|---|
| [Diagonal blades](blades-1161974.json) | Nearest linked tube | 1 | 26 | $\log(2)/26$ |
| [Four curled fins](fins-394727.json) | Single tube | 2 | 104 | $\log(2)/104$ |
| [Helical chain](helix-2189166.json) | Blended linked tubes | 5 | 295 | $\log(2)/295$ |
| [Stepped arches](arches-303603.json) | Single tube | 3 | 321 | $\log(2)/321$ |

The [two-evaluator audit](audit.json) also compares these periodic points and intermediate numerical excursion states with the 250,000-point forward portraits. For Four curled fins, every sampled excursion state lies within 0.018 RMS radii of the portrait. The 99th percentile of nearest distances from the later half of the same portrait to its earlier half is 0.052 RMS radii.

The following argument details the helical-chain packet. The other three use their stored primitive period, chart radii and block lengths in the same inequalities. For Case 2189166 there is a compact invariant set on which the topological entropy of the original map is at least $`\log(2)/295`$. Both interval calculations verify the same two inverse branches of the 295th iterate. The measured helical point cloud comes from a separate forward orbit. Its equality with this subsystem and its asymptotic attraction remain open. [Model and parameter set](../case-2189166.json) · [Certificate](helix-2189166.json)

## The periodic point and its local inverse

Write $`F`$ for the blended map and $`H=F^5`$. An interval Krawczyk inclusion of radius $`10^{-115}`$ encloses a unique zero $`x_*`$ of $`H(x)-x`$. Disjoint enclosures exclude period one. All three multipliers are outside the unit disk by the interval Routh test described in the [study report](../README.md). Coordinates below use the exact stored real basis $`B`$ and its rational inverse.

Use the infinity norm in the coordinates $`u=B^{-1}(x-x_*)`$. The closed box $`W=\{\|u\|_\infty\le10^{-15}\}`$ lies in a smooth domain for all five updates. For the chart return $`h`$ the interval Jacobian enclosure and a stored preconditioner $`C`$ give

```math
 q=\sup_W\|I-C\,Dh\|_\infty,\qquad
 q+\|C\|_\infty<1,\qquad
 k=\frac{\|C\|_\infty}{1-q}<0.703.
```

For each target $`v\in W`$ the function $`u\mapsto u-C(h(u)-v)`$ maps the closed box into itself and contracts by at most $`q`$. The fixed point defines a unique inverse branch $`G:W\to W`$ with $`G(0)=0`$. Subtracting the fixed-point equations for two targets yields the displayed Lipschitz bound $`k`$. The verified Jacobian condition also implies invertibility of $`C`$ and of every covered derivative.

## A second branch on the same target

Let $`V`$ be the concentric chart box of radius $`5\times10^{-18}`$. The return itinerary consists of 59 regular inverse blocks of length five. In each block, a parameterized Krawczyk inclusion maps every point of the preceding target enclosure to a unique point in the new source enclosure. The first target is a chart box of radius $`10^{-17}`$ around the stored decimal center. It contains $`V`$ for every possible location of $`x_*`$ inside its interval enclosure. Composing the 59 continuous local inverses defines a second branch $`J`$ of $`H^{59}=F^{295}`$ on $`V`$.

The largest Krawczyk inclusion ratio is below 0.975 with either evaluator. The image $`J(V)`$ stays within $`4.590\times10^{-18}`$ of the periodic point in chart units. At least one coordinate interval lies strictly outside the box of radius $`k^{59}(5\times10^{-18})`$. Thus $`J(V)\subset V`$ and $`J(V)`$ is disjoint from $`G^{59}(V)`$. All inequalities in the checker use outward intervals, including the image separation. Decimal values here are rounded for display.

## Entropy

The two continuous inverse branches map the compact target $`V`$ into disjoint compact subsets of itself. The unions of all branch images of word length $`n`$ form a nested sequence of nonempty compact sets. Their intersection $`K`$ is invariant under $`F^{295}`$. Every finite binary itinerary has an extension in $`K`$. Points with different itineraries separate by at least the positive distance between the two first-level images at the first different symbol. The entropy of $`F^{295}|_K`$ is at least $`\log 2`$. Taking the union of the first 295 forward images of $`K`$ gives a compact $`F`$-invariant set and divides the bound by 295. The validated regular domains cover all intermediate images in this argument.

## Full two-symbol shift and expansion

Both inverse branches of the return are contractions in the same chart norm. For a parameterized inverse block with preconditioner $C_i$ and verified residual row norm $q_i<1$, subtracting the two fixed-point equations gives a Lipschitz constant at most $\lVert C_i\rVert_\infty/(1-q_i)$. Multiplying these outward bounds along the excursion gives a strict contraction. The other branch has bound $k^{m}$ for the corresponding number of local inverse returns. The checker explicitly requires both bounds below one. [Base-case bounds from both evaluators](audit.json) · [Twelve follow-up cases](../three/snapback/audit.json)

The nested target sets for any infinite binary word therefore have diameters tending uniformly to zero. Each word gives one point, distinct words give distinct points by the disjoint first-level images, and the forward return deletes the first symbol. This coding is a homeomorphism from the one-sided binary shift to its compact invariant Cantor set. The bound on inverse derivatives also gives uniform expansion for the return along that set in the chart norm. In particular, every periodic binary itinerary is realized. For Diagonal blades the larger of the two branch bounds is below $2.59\times10^{-7}$ in both implementations. [Diagonal blades packet](blades-1161974.json)

This construction uses regular returns to a repeller. The explicit local inverse test addresses the expanding-neighborhood condition in the corrected snap-back-repeller literature. [Marotto, 2005](https://www.sciencedirect.com/science/article/pii/S0960077904006563)

## Recheck

The packet encodes the exact binary input parameters through their stored decimal binary64 representations, rational chart matrices, decimal enclosure centers, preconditioners and radii. The two evaluators share the certificate, mpmath interval arithmetic and Krawczyk driver. The first computes derivatives by forward automatic differentiation of nested reflectors. The second uses closed-form Gaussian posterior covariance and both principal torus tangent directions. Their upper bounds on the local inverse Lipschitz constant are 0.7020 and 0.6731 respectively. Independent code paths share the numerical arithmetic library. [Checker](../../../../research/shape_search_20261121/projection/snapback.py) · [Analytic interval implementation](../../../../research/shape_search_20261121/projection/analytic_bounds.py)

Run from the repository root.

```powershell
$src = 'research/shape_search_20261121/projection'
$res = 'results/shape_search_20261121/projection'
python "$src/snapback.py" --verify --engine jet --record "$res/case-2189166.json" --out "$res/snapback/helix-2189166.json"
python "$src/snapback.py" --verify --engine analytic --record "$res/case-2189166.json" --out "$res/snapback/helix-2189166.json"
```

Run the full audit for all four packets.

```powershell
python "$src/audit_entropy.py" --source $res
```

The numerical proposals record how the chains were initialized. The verifier reads the frozen interval certificates directly. [Diagonal blades proposal](proposal-1161974.json) · [Four curled fins proposal](proposal-394727.json) · [Helical chain proposal](proposal-2189166.json) · [Stepped arches proposal](proposal-303603.json)

The inequalities are strict, involve finitely many boxes, and stay on continuous differentiable branches under sufficiently small parameter perturbations. The implicit function theorem continues each hyperbolic periodic root. These facts give an open parameter neighborhood with the same symbolic lower bound. This continuity argument does not give the size of that neighborhood.

Selected packets also verify closed intervals of the original relaxation parameter. Holding all other binary parameters fixed, the two interval calculations enclose a continued repelling root and both branches for every real relaxation within the stored absolute radius of the center. The explicit guaranteed radii range here from $10^{-78}$ to $10^{-18}$. These small sufficient radii arise from using frozen branch centers throughout the return. [Parameter bands](../bands/) · [Uniform parameter checker](../../../../research/shape_search_20261121/projection/entropy_band.py)

# Pinched graphs over a two-torus

The 1,536 cases vary two or three integer wave vectors, their phase offsets and the multiplier of a quasiperiodically forced map. The search selected 36 geometries by farthest-point sampling of distance and height quantiles. The dense sample of each case records 240,000 iterates after an 8,000-step transient. Three initial conditions independently checked the typical fiber rate. The analytic results below apply to every case, including the ones omitted by the visual selection. Pinched skew products and strange nonchaotic attractors have an established literature. [Parameters and numerical receipts](../../../results/shape_search_20261121/pinched/manifest.json) · [Jäger, *On the structure of strange non-chaotic attractors in pinched skew products*](https://arxiv.org/abs/math/0402360)

## Dynamics and proof

The base point $`q=(\theta,\phi)`$ advances by $`Tq=q+(\sqrt2-1,\sqrt3-1)`$ modulo one. Integer wave vectors are nonzero and primitive. The height map preserves the unit interval.

```math
a(q)=\kappa\prod_{j=1}^m\sin^2\!\bigl(\pi(k_j\theta+l_j\phi+b_j/32)\bigr),\qquad
h'=\frac{a(q)h}{1+(\kappa-1)h}.
```

The uniquely ergodic base translation preserves Haar probability. The integrable function $`\log a`$ has mean $`L=\log(\kappa/4^m)>0`$ because every integer phase projects Haar measure to uniform measure on the circle. The search takes $`\kappa/4^m`$ from 1.5, 2, 3 and 4. The zero graph has fiber exponent $`L`$ at typical base points.

In the reciprocal coordinate $`y=1/h`$, a positive orbit satisfies $`y'=(y+\kappa-1)/a(q)`$. Pulling back the top of the interval gives the upper invariant graph $`h_*=1/y_*`$ with

```math
y_*(q)=(\kappa-1)\sum_{j=1}^\infty
      \left(\prod_{i=1}^j a(T^{-i}q)\right)^{-1}
```

at almost every base point. The ergodic theorem gives exponential decay of the summands and of the omitted initial-condition term. Everywhere, the pullback of the top defines an upper semicontinuous graph since the height maps are increasing and the pullback sequence decreases. At a pinch the next height is exactly zero. All forward translates of each pinching curve therefore have zero graph height, and their union is dense. The upper graph is discontinuous everywhere it is positive, a set of full base measure.

On typical positive orbits, differences of reciprocal heights contract at exponential rate $`-L`$. The heights lie in $`[0,1]`$, so this also gives forward attraction in the height coordinate for positive starts. For the exponent on the upper graph, its logarithmic increment is $`\log a-\log(1+(\kappa-1)h_*)`$, an integrable coboundary. Its integral vanishes by telescoping and recurrence to a positive-measure set where $`\log h_*`$ is bounded. Averaging the logarithm of the fiber derivative consequently gives $`-L`$. The base contributes two zero exponents. The typical spectrum is thus $`(0,0,-L)`$. These statements concern the measurable attracting graph and its typical dynamics.

## What the portrait measures

The two base angles are drawn as toroidal and poloidal angles. Major radius is 1.32, and minor radius is $`0.20+0.95h`$. The lower graph appears as the surface of a thin tube. Large excursions in height form ridges in the plotted sample. These coordinates describe a drawing of the full state, with a chosen tube-radius scale. Topological equivalence between different parameter choices and dimensions of graph closures were outside the calculations.

The forward solver evolves log height. Three of the 25,920,000 displayed heights rounded to zero after conversion back to floating point. The log-height orbit and the fiber derivative continued through these events. Across 108 long runs the largest discrepancy from the predicted fiber exponent was below 0.00684 per iteration.

An independent 180-bit reciprocal-series calculation checked two prescribed base phases per case at pullback depths 256, 1,024 and 2,048. At depth 2,048, the largest discrepancy with forward computation in log height was below 1.96 × 10⁻⁸ over all 72 comparisons. [Numerical pullback checks](../../../results/shape_search_20261121/pinched/pullback-check.json)

## Reproduction

NumPy and Numba run the search and sampling. A single case can be replayed directly from the saved parameter index. RUN is an external scratch directory.

```text
python research/shape_search_20261121/pinched/pullback.py results/shape_search_20261121/pinched/manifest.json --out RUN/pullbacks.json
python research/shape_search_20261121/pinched/search.py --count 1536 --short 16000 --out RUN/screen
python research/shape_search_20261121/pinched/search.py --replay results/shape_search_20261121/pinched/manifest.json --cases 156 1267 --steps 240000 --out RUN/clouds
```

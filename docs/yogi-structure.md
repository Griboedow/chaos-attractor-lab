# Elementary structure of the scalar-quadratic Yogi map

This note accompanies the [numerical study](yogi-geometry.md). The algorithm is due to [Zaheer et al. (2018)](https://www.satyenkale.com/papers/yogi.pdf). The calculations below explain implementation controls and constrain possible interpretations. They are not claimed as historically new theorems or proofs of chaos.

Write `0<b,c<1`, `e>0`, and

```text
s = sign(v-x^2),  M=b*m+(1-b)*x,
V=v-(1-c)*s*x^2,  X=x-M/(sqrt(V)+e).
```

## Positivity, reflection, and volume

For `v>=0`, `V>=0`. More specifically,

```text
V >= c*v,
V >= delta*x^2,  delta=min(c,1-c).
```

The first inequality is already used in the original paper's Appendix B. For the second, the branch `v>x^2` gives `V>=c*x^2`, the branch `v<x^2` gives `V>=(1-c)*x^2`, and equality gives `V=x^2`.

The reflection `R(x,m,v)=(-x,-m,v)` commutes with the update. Reflected lobes or cold-start pairs must therefore not be counted as unrelated discoveries.

Away from the switching surface and `V=0`, first replace `(m,v)` by `(M,V)` while holding `x` fixed. This triangular transformation has determinant `b`. Then replace `x` by `x-M/(sqrt(V)+e)` while holding `M,V` fixed; that shear has determinant1. Thus `det DF=b` on every smooth branch. This is a local volume contraction, not a global attraction or boundedness theorem.

## A genuine jump and overlapping inverse branches

At `v=x^2`, `x!=0`, the upper-side and lower-side limits of `V` are respectively `c*x^2` and `(2-c)*x^2`, whereas the exact `sign(0)=0` value is `x^2`. The one-sided `X` values are

```text
x-M/(sqrt(c)*abs(x)+e),
x-M/(sqrt(2-c)*abs(x)+e).
```

Even when these two `X` values happen to agree, the third coordinate jumps by `2*(1-c)*x^2`. A global homeomorphic conjugacy with a continuous map on the full state domain would preserve continuity and is therefore impossible. This observation does **not** rule out equivalence of selected invariant subsets.

Given a target `(X,M,V)`, recover

```text
x = X+M/(sqrt(V)+e),
m = (M-(1-b)*x)/b,
v_s = V+(1-c)*s*x^2.
```

For `x!=0`, both strict branches `s=±1` are valid, distinct, nonnegative predecessors whenever

```text
max(c,1-c)*x^2 < V < (2-c)*x^2.
```

Hence the full map is noninjective, despite its nonzero branch determinant. The frozen controls reconstruct 512 such branch inverses with maximum residual below `9e-16`. Global noninjectivity does not establish noninjectivity after restriction to an observed attractor.

## The fixed line and competing convergence

Every `(0,0,v)` with `v>=0` is fixed. For `v>0`, put `d=sqrt(v)+e`. The transverse linear block is

```text
[ 1-(1-b)/d    -b/d ]
[   1-b          b  ]
```

Its determinant is `b`, its trace is `1+b-(1-b)/d`, and the third multiplier is1. The quadratic Schur conditions reduce to

```text
d > (1-b)/(2*(1+b)).
```

Thus sufficiently large memory states are linearly stable in the two transverse directions, with a neutral direction along the line. This explains the possibility of competing fixed outcomes; the linear calculation alone does not prove their basins. It also explains why a zero leading rate at a settled state need not mean quasiperiodicity.

## A bound for a restricted zero-momentum initialization

Assume `m0=0`, `v0>=0`, `x0` finite, and `b<sqrt(c)`. Define

```text
delta = min(c,1-c),  r=b/sqrt(c),
K = (1-b)/(sqrt(delta)*(1-r)),
D = b*K/(1-b),
B = max(abs(x0),D+K).
```

Then every autonomous iterate satisfies

```text
abs(x_t)<=B,  abs(m_t)<=B,
0<=v_t<=max(v0,(2-c)*B^2).
```

The same bound holds if both moments are bias-corrected on every update, still with **zero initial first moment**. It does not cover the installed Optax release's nonzero initial first moment, nor the high-`b` sunburst and two-direction cases. It covers the paired-fan parameters.

To prove it, write `m_{t+1}=(1-b)*sum_{k=0}^t b^k*x_{t-k}`. Positivity gives `v_{t+1}>=delta*c^k*x_{t-k}^2`. The resulting geometric-series estimate bounds each autonomous displacement by `K`; if the denominator vanishes before adding `e`, all relevant `x` terms vanish and the step is zero. With bias correction the corresponding finite sum is multiplied by

```text
sqrt(1-c^(t+1))*(1-r^(t+1))/(1-b^(t+1)) <= 1,
```

so the same step bound applies. The normalized momentum `mhat=m_{t+1}/(1-b^(t+1))` is a convex average of recent `x` values. Their separation from `x_t` is at most `k*K`; the truncated geometric distribution has mean age at most `b/(1-b)`. Therefore `abs(mhat-x_t)<=D`.

For `x_t>D`, the update points left; for `x_t<-D`, it points right. If `abs(x_t)<=D`, a step can reach at most `D+K`. If an inward step crosses zero from outside this interval, its overshoot is also bounded by `K`. Induction gives the bound for `x`. The raw momentum is a weighted sum with total weight at most1, giving its bound. Finally, when `v_t>x_t^2` the second memory decreases; otherwise its next value is at most `(2-c)*B^2`.

This is a scalar specialization of familiar momentum/accumulator estimates, not a priority claim. Related Adam bounds under `beta2>beta1^2` appear in [Dereich et al. (2025)](https://arxiv.org/abs/2509.10476). Their preconditioner hypotheses must be checked before applying those results to Yogi; the recurrences are not interchangeable. Our 512 finite numerical bound controls test the implementation, while the argument above supplies the mathematical statement.

## What remains unproved

None of these identities proves a compact trapping region or chaotic attraction for cases 7590 or4214. They do not identify a Markov partition, a uniformly expanding bundle, a physical invariant measure, or a new topological class. In particular, two positive finite branch rates remain a numerical observation, with independent finite-separation checks as corroboration. See the [frozen structural controls](../results/yogi-geometry/structure/analysis.json) and [priority audit](yogi-priority.md).

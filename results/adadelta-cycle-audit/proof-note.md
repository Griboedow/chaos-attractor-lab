# AdaDelta: an apparent two-cycle and the effect of positive epsilon

A close prior-art check reveals a material distinction between a zero-epsilon family of exact cycles and a positive-epsilon trajectory that slowly drifts. For the scalar quadratic AdaDelta recurrence printed below, a nontrivial period-two orbit with a positive denominator regularizer exists **if and only if `alpha*curvature=2`**. At that critical product it is a family of symmetric cycles. This is an exact algebraic statement about the full three-variable map, not a numerical inference from a long trace.

The check concerns [Sebastian Bock's dissertation, DOI 10.5283/epub.58326](https://epub.uni-regensburg.de/58326/), equation 5.2 and Section 5.2.3, including Theorem 5.4 and Figure 52. Typeset printed pages 91 and 101-102 were inspected. The displayed recurrence and the source's nonzero-epsilon example were read directly. This audit is narrow: it does not provide a global convergence theorem for AdaDelta or change the separate evidence for Adam and Yogi chaos.

## Exact recurrence and conventions

Let the objective be `f(x)=k*x^2/2`, with `k>0`. Let `alpha>0`, `0<=beta<1`, and `E>0`. The source writes its regularizer as `epsilon^2`; we write `E` to distinguish it from the epsilon conventions used by other optimizers. Let `v>=0` be the squared-gradient memory and `d>=0` the squared-update memory. One update is

$$
V=\beta v+(1-\beta)k^2x^2,\qquad
\delta=kx\sqrt{\frac{d+E}{V+E}},
$$

$$
D=\beta d+(1-\beta)\delta^2,\qquad
X=x-\alpha\delta.
$$

The update-memory recurrence uses `delta^2`, before multiplication by `alpha`. Replacing it with `(alpha*delta)^2` would define a different algorithm and invalidate the following parameter condition. The inspected actual PyTorch AdaDelta implementation uses the stated convention, with its `eps` equal to `E`.

## Classification of all nontrivial period-two orbits for positive E

Suppose a full-state period-two orbit exists. If its two positions are equal, the position update forces both to be zero; the two memory recurrences then force both memories to zero. This is the fixed point, not a nontrivial two-cycle.

Otherwise the two positions must have opposite signs. An update from a positive position decreases it, and an update from a negative position increases it, because the square-root ratio is positive. Write the positions as `a` and `-b`, where `a,b>0`.

Both unscaled updates have the same magnitude, `(a+b)/alpha`. The two-step update-memory recurrence consequently forces

$$
d_0=d_1=\frac{(a+b)^2}{\alpha^2}.
$$

The squared-gradient memories used in the two updates are

$$
V_a=\frac{k^2(a^2+\beta b^2)}{1+\beta},\qquad
V_b=\frac{k^2(\beta a^2+b^2)}{1+\beta}.
$$

Squaring the two position-update equations and cancelling their common positive factors gives

$$
(V_a+E)b^2=(V_b+E)a^2.
$$

Substitution factors this equation as

$$
(b^2-a^2)\left[\frac{\beta k^2(a^2+b^2)}{1+\beta}+E\right]=0.
$$

The bracket is strictly positive. Thus `a=b`: symmetry is a consequence, not an assumed restriction on the possible cycles. The two memories reduce to

$$
v=k^2a^2,\qquad d=4a^2/\alpha^2.
$$

The remaining position equation becomes

$$
4(k^2a^2+E)=\alpha^2k^2(4a^2/\alpha^2+E),
$$

or `E*(4-alpha^2*k^2)=0`. Since `E>0`, `alpha>0` and `k>0`, this proves necessity of `alpha*k=2`.

Conversely, if `alpha*k=2`, every `a>0` gives the exact orbit

$$
(a,k^2a^2,k^2a^2)\;\longleftrightarrow\;
(-a,k^2a^2,k^2a^2),
$$

where the coordinates are `(x,d,v)`. Direct substitution proves sufficiency. This completes the classification for the stated positive-regularizer recurrence. When `E=0`, the familiar symmetric family with `d=4a^2/alpha^2` exists for arbitrary positive `alpha,k`; the positive-`E` conclusion must not be transferred to that singular case.

## Replaying the source's example

The displayed example uses `alpha=0.1`, `beta=0.9`, `k=1`, source `epsilon=1e-4` and initial `(x,d,v)=(1,400,1)`. Thus `E=1e-8` and `alpha*k=0.1`, outside the exact period-two condition.

After two ordinary binary64 updates, the displacement from the initial full state is approximately

```text
(-1.89525e-8, -1.51620e-6, -1.99500e-9).
```

Independent 70- and 90-digit evaluations using the exact rational parameters agree that, after 5,000 updates,

```text
x = 0.999975037252175489...
d = 399.980067955994009...
v = 0.999950174877484781...
```

The actual PyTorch 2.15.0 CPU optimizer, with these memories explicitly initialized and `eps=1e-8`, agrees at that endpoint to a maximum state discrepancy of `4.84e-12`. Its complete 5,000-step output is retained. A source-like short plot cannot visually resolve the small amplitude change.

A separate binary64 continuation to 100,000,000 updates reaches `x≈0.05000248`, `d≈1.00013717`, `v≈0.00250035`. The oscillation amplitude has substantially changed. Controls retain the exact zero-`E` cycle, retain the positive-`E` critical-product cycle, show outward finite drift at `alpha=3`, and reach the floating-point underflow scale with the larger `E=1e-6` at `alpha=0.1`. That last numerical outcome is not presented as an exact-real stopping or convergence proof. Every declared control is retained.

## Consequence for the literature comparison

The source's zero-epsilon cycle family is reproduced. Its continuation argument treats the squared-gradient memory as a fixed parameter in a reduced system. A root of that reduced system must still satisfy the omitted memory return equation to be a cycle of the full map. The classification above shows that a nontrivial positive-`E` cycle cannot satisfy all three return equations away from `alpha*k=2`.

Accordingly, this repository credits the zero-epsilon family as prior work but does not rely on the quoted positive-epsilon extension as evidence of an exact full-state cycle. The stated source example is a slowly changing oscillation in the reproduced recurrence. This correction concerns that claim and convention; absence of a two-cycle does not rule out other forms of nonconvergence.

The audit was motivated by the search for unusual attractors: a nearly closed portrait and a long flat-looking amplitude trace are insufficient. The [Yogi report](yogi-parameter-plane.md) keeps its separate long-run, precision, library and initial-state checks. No new attractor is claimed from the AdaDelta replay.

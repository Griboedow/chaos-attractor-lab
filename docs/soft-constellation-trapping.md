# A trapping certificate for the Gaussian-posterior folded band

**Status, 2026-11-01:** an outward interval certificate establishes a compact
attracting invariant set for Case 156 of the
[Gaussian-posterior circle experiment](soft-constellations.md). The original
initial state, and a small explicit neighborhood of it, enter the trapping
region after eleven iterations. Chaotic dynamics on the limiting set and
historical uniqueness are not established by this certificate.

The [published portrait and rotation](soft-constellations.md) depict numerical
states, rather than a rendering of a proved complete invariant set. The
certificate supplies an attracting region containing the entire terminal
100,000-state reference sample. Equality of that numerical ribbon with the
maximal invariant set below is not asserted.

## Exact model and region

The model uses six fixed sites, the unit circle in the xy-plane, and

```math
p_i(x)=\frac{\exp(-\|x-b_i\|^2/\tau)}
{\sum_j\exp(-\|x-b_j\|^2/\tau)},\qquad m(x)=\sum_i p_i(x)b_i,
```

```math
v(x)=2m(x)-x,\qquad r(x)=\sqrt{v_1(x)^2+v_2(x)^2},\qquad
F(x)=\beta x+\beta\frac{(v_1(x),v_2(x),0)}{r(x)}+(1-2\beta)m(x).
```

Here `beta=0.970660805775784` and `tau=0.03731678183231506` denote the
**exact real lifts of the stored binary64 values**, as do all site
coordinates. These are the unchanged parameters of Case 156; idealized
polyhedral symmetry or rounded decimal replacements are not used.

The complete cubical region and parameters are retained in
[proposal.npz](../results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz),
SHA-256 `1ae5e76bde96c9193eaca526eb6d233ddd5a68740311e4391ae6944d661957eb`.
For its stored origin `o`, shape `(n_x,n_y,n_z)` and occupied addresses `S`,
decode `k=(i*n_y+j)*n_z+l` and define

```math
Q=\bigcup_{(i,j,l)\in S}
\prod_{a=1}^3[o_a+h s_a,\ o_a+h(s_a+1)],\qquad h=2^{-8}.
```

The cubes are closed and their coordinates are exact dyadic rationals.
There are **35,784 cubes**. Their union lies within

```math
[0.33203125,0.52734375]\times[-1.06640625,-0.6171875]
\times[-0.359375,0.21875].
```

The certificate checks the union, not the larger rectangular hull.

## Certified inclusion and its consequence

The complete computation establishes

```math
\boxed{F(Q)+[-2^{-35},2^{-35}]^3\subset Q},\qquad
\boxed{r(x)>0.6697\quad(x\in Q)}.
```

The radius bound excludes the projection singularity throughout the region.
Consequently, `F` is smooth on a neighborhood of `Q`, and the positive
coordinate margin implies `F(Q) subset int(Q)`.

**Proposition.** For this exact model and region,

```math
K=\bigcap_{n\ge0}F^n(Q)
```

is nonempty, compact, satisfies `F(K)=K`, and attracts a neighborhood of `Q`.

**Proof.** The nonempty compact sets `F^n(Q)` are nested. Their intersection
is therefore nonempty and compact. Continuity and compactness imply that
the image of this nested intersection is the intersection of its images,
so `F(K)=K`. Nested compact sets approach their intersection in Hausdorff
distance. Thus `sup_{x in Q} dist(F^n(x),K)` tends to zero. Since the image
of `Q` lies in its interior, continuity supplies an open neighborhood of
`Q` whose image lies in `Q`; this neighborhood is attracted to `K` as well.

This conclusion neither selects a unique invariant probability measure nor
establishes transitivity, positive entropy, or a positive asymptotic
Lyapunov exponent on `K`. It does not exclude periodic or other invariant
subsets within `K`.

## How the inclusion is checked

Each source cube is either accepted whole or subdivided into eight closed
dyadic subcubes. For a subcube `X`, center `c`, and chosen binary64 matrix
`A`, the mean-value formula gives the enclosure

```math
F(X)\subset F(c)+A(X-c)+(DF(X)-A)(X-c).
```

The code encloses `F(c)` and `DF(X)` with outward interval arithmetic.
It adds the exact margin `2^-35` in every output coordinate. An outward
separating-axis calculation rejects target cubes only when they cannot
intersect this enclosure. Every remaining target cube must belong to `Q`.
Failure to resolve a subcube is a failed certificate check, never evidence
of an escaping orbit.

The Jacobian extension uses `Dm=(2/tau) Cov_p(b)` and the exact formula in
the [model report](soft-constellations.md#model-and-known-ingredients).
Posterior weights are evaluated with pairwise logit differences relative
to one site. All basic operations are expanded with `nextafter`; no
fast-math transformations or parallel floating-point reductions are used.
The square-root enclosure has the same binary64 arithmetic assumptions
as the earlier [circle trapping certificate](raar-trapping/manuscript.md).

The exponential requires no error assumption for a library exponential.
Arguments are reduced by exact powers of two to `[-1/8,1/8]`. A degree-12
Taylor polynomial is evaluated with outward arithmetic, with the absolute
remainder bounded by

```math
e^{1/8}\frac{(1/8)^{13}}{13!}
\le\frac87\frac{(1/8)^{13}}{13!}<2^{-71}.
```

Repeated interval squaring reverses the range reduction. Monotonicity gives
the enclosure for an interval argument from its endpoint enclosures.

The completed audit covers **1,489,446 accepted subcubes**, with **207,666
subdivisions** and maximum depth **4**. The exact accounting identity is
`1,489,446 = 35,784 + 7*207,666`. Every root cube is processed once;
there are no unresolved leaves. Reversing the root order and changing the
batch size from 1,000 to 750 reproduces the aggregate counts, depth and
minimum radius bound exactly.

Independent 100-digit tests check 1,031 exponential arguments, 384 state
evaluations and 3,456 derivative entries against a separately expressed
nested-reflector formula. These tests support the implementation; the
inclusion itself rests on exhaustive interval coverage of `Q`.

## Entry of the original initial state

Let `x0=(0.317,-0.211,0.419)`, again interpreted as exact binary64 lifts.
An interval propagation proves

```math
F^{11}(x_0+[-2^{-35},2^{-35}]^3)\subset Q.
```

Thus this explicit neighborhood of the original start approaches `K`.
Natural and mean-value interval enclosures are intersected at each step.
The final enclosure has the following recorded binary64 endpoints
(their exact bits, rather than rounded decimal real replacements, are used):

```text
x: [ 0.37964081635036795,  0.37964081690007070]
y: [-0.85444580926089770, -0.85444580771712990]
z: [-0.14198734332327545, -0.14198734257609386]
```

An independent 80-digit mpmath interval calculation uses the nested
reflection formula and confirms eleven-step entry for the smaller initial
radius `2^-160`. Its natural interval propagation at radius `2^-35` became
unbounded through overestimation; that failed control is retained. It does
not contradict the tighter mean-value certificate for the larger box.

## Selection, failed attempts and interpretation

The region was proposed from the frozen numerical sample using sampled
corner remainders. Those heuristic bounds are not used as proof. The
completed interval audit checks the resulting fixed union independently
of the proposal's claimed closure.

Earlier Case 29 proposals at spacings `1/128`, `1/256`, `1/512` and `1/1024`
did not close: two failed the proposed alternation enclosure, one left its
declared bounding domain, and one exceeded the new-cell budget. These
attempts and an initial implementation failure are retained. None establishes
nonexistence of a trapping region for Case 29.

The known Gaussian mean-shift/MMSE/proximal ingredients and earlier
reflection dynamics are documented in the
[priority review](soft-constellation-priority.md). Two-band alternation
is also established prior art. The present additional result is a
reproducible trapping and entry certificate for this parameter set;
historical priority and a proof of chaotic attraction remain open.

[Evidence packet and reproduction](../results/soft-trapping/README.md).

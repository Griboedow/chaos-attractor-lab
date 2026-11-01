# An expanding-circle factor in cubic FastICA

The standard one-unit cubic FastICA update admits an explicit finite-data
construction with a compact attracting set of topological entropy at least
`log(3)` per iteration. The statement holds for every `0 <= epsilon <= 1/25`.
It follows from exact moment identities, a strict trapping tube and an explicit
continuous factor map onto angle tripling. It does not depend on a finite-time
Lyapunov estimate.

The observations are deliberately constructed and whitened; an identifiable
independent-source ICA model is not assumed. The algorithm and the expanding
circle mechanism are established mathematics. Historical priority of this
particular finite-data realization remains unresolved.

[Main-page portraits](../README.md#structured-data-fastica) ·
[Evidence and reproduction](../results/fastica-tripling/README.md) ·
[Earlier empirical FastICA study](fastica-dynamics.md)

## Data and update

Use the following 20 equally weighted base observations `(x,y)`: four copies of
each of `(1,0), (-1,0), (0,1), (0,-1)`, and one copy of each of
`(sqrt(6),0), (-sqrt(6),0), (0,sqrt(6)), (0,-sqrt(6))`.
Set `R=x^2+y^2`, `a=x(R-4)`, `b=y(R-4)` and
`sigma=sqrt(1-6 epsilon^2)`. For each base observation, take the Cartesian
product of two six-element lists
`[sqrt(3), -sqrt(3), 0, 0, 0, 0]`, denoting their entries by `xi, eta`.
The resulting **720 equally weighted rows** are

$$
X=(x,y,\epsilon a+\sigma\xi,\epsilon b+\sigma\eta).
$$

The data are fixed throughout iteration. In exact real arithmetic,
`E[X]=0` and `E[XX^T]=I`. The construction is defined for
`0 <= epsilon < 1/sqrt(6)`; the theorem below uses the smaller interval.
There is no time-dependent forcing, resampling, added memory or modified update.

For a unit weight `w=(A,B,C,D)` in `S^3`, cubic FastICA uses

$$
h(w)=\mathbb E[X(X\cdot w)^3]-3\mathbb E[(X\cdot w)^2]w,
\qquad F(w)=\frac{h(w)}{\|h(w)\|}.
$$

With exact whitening, this equals `h=grad(P)/4`, where the homogeneous
fourth-cumulant polynomial is

$$
\begin{aligned}
P={}&A^4+B^4-6A^2B^2+24\epsilon(A^3C+B^3D)\\
 &+\epsilon^2[72(A^2C^2+B^2D^2)-36(A^2D^2+B^2C^2)]\\
 &+72\epsilon^3(AC^3+BD^3)
 -18\epsilon^4(C^4+D^4)-216\epsilon^4C^2D^2.
\end{aligned}
$$

For an independent derivation, put
`V=x[A+epsilon C(R-4)]+y[B+epsilon D(R-4)]`. The base moments give

$$
\begin{aligned}
\mathbb E[V^2]&=A^2+B^2+6\epsilon^2(C^2+D^2),\\
\mathbb E[V^4]&=\frac25[(A-3\epsilon C)^4+(B-3\epsilon D)^4]
 +\frac{18}5[(A+2\epsilon C)^4+(B+2\epsilon D)^4].
\end{aligned}
$$

Each noise list has moments `(1,0,1,0,3)` through order four. Thus, with
`N2=(1-6 epsilon^2)(C^2+D^2)`,
`P=E[V^4]+6 E[V^2] N2+3 N2^2-3 ||w||^4`.
[The symbolic check](../results/fastica-tripling/code/exact_check.py) verifies
these identities and the rational inequalities below without floating-point
rounding.

## Attracting set and entropy

**Proposition.** For every `0 <= epsilon <= 1/25`, let

$$
T=\{(A,B,C,D)\in S^3:\sqrt{C^2+D^2}\leq \tfrac25\sqrt{A^2+B^2}\}.
$$

Then `F` is defined on `T`, `F(T)` lies strictly inside `T`, and

$$
K=\bigcap_{n\geq0}F^n(T)
$$

is nonempty, compact and invariant, and attracts `T` uniformly in distance.
The restricted map `F|K` has a continuous onto factor
`theta -> -3 theta (mod 2 pi)`. Consequently,
`h_top(F|K) >= log(3)`. The same lower bound holds after identifying
the equivalent ICA weights `w` and `-w`.

### Strict trapping

Write `r=sqrt(A^2+B^2)` and `s=sqrt(C^2+D^2)`. The first two components at
`epsilon=0` are

$$
G=(A^3-3AB^2,\ B^3-3BA^2),\qquad \|G\|=r^3.
$$

In complex notation `G_1+iG_2=conjugate(A+iB)^3`. Direct differentiation
of `P` and the triangle inequality give

$$
\begin{aligned}
\|h_{AB}-G\|&\leq18\epsilon r^2s+36\epsilon^2rs^2+18\epsilon^3s^3,\\
\|h_{CD}\|&\leq6\epsilon r^3+36\epsilon^2r^2s+54\epsilon^3rs^2
 +108\epsilon^4s^3.
\end{aligned}
$$

For example, the quadratic-in-`epsilon` terms in `h_AB` are
`(36AC^2-18AD^2, 36BD^2-18BC^2)`, whose norm is at most `36 r s^2`.
The quartic terms in `h_CD` are
`-18(C(C^2+6D^2), D(D^2+6C^2))`, bounded by `108 s^3`.
The other terms have the displayed bounds by the same coordinate norm
inequalities.

For `s/r <= k=2/5` define

$$
\begin{aligned}
L(\epsilon)&=1-18\epsilon k-36\epsilon^2k^2-18\epsilon^3k^3,\\
U(\epsilon)&=6\epsilon+36\epsilon^2k+54\epsilon^3k^2+108\epsilon^4k^3.
\end{aligned}
$$

`L` decreases and `U` increases on the parameter interval. At its endpoint,

$$
L_* = \frac{1372481}{1953125},\qquad
U_* = \frac{12871614}{48828125},\qquad
kL_*-U_* = \frac{853196}{48828125}>0.
$$

Hence `||h_AB|| >= r^3 L_* > 0` and

$$
\frac{\|F_{CD}\|}{\|F_{AB}\|}
\leq\frac{U_*}{L_*}
=\frac{12871614}{34312025}
\simeq0.3751342<0.4.
$$

The normalization cancels in this ratio. This proves strict trapping and
excludes a zero numerator everywhere in `T`. The nested nonempty compact
sets `F^n(T)` have intersection `K`; continuity gives `F(K)=K`.
Compactness also implies `sup_{z in F^n(T)} dist(z,K) -> 0`.
Since `K` lies in the interior of `T`, this is attraction of a neighborhood.

### Explicit circle factor

Identify `T` with `S^1 x closed_disk(2/5)` by
`w=(cos(theta),sin(theta),u_1,u_2)/sqrt(1+||u||^2)`.
The relative complex number `(h_A+i h_B)/(G_1+i G_2)` has positive real
part because its distance from 1 is at most `1-L_*<1`. Its argument
`q(w)` therefore has a continuous real branch in `(-pi/2,pi/2)`.
The angular component of `F` satisfies

$$
\theta(Fw)=-3\theta(w)+q(w)\pmod{2\pi}.
$$

The uniformly convergent series

$$
H(w)=\theta(w)+\sum_{j=0}^{\infty}(-3)^{-j-1}q(F^jw)\pmod{2\pi}
$$

defines a continuous map with `H(Fw)=-3H(w)`.
For each fixed `u`, the correction is periodic in `theta`, so `H` has
degree one on that circle and is onto. Thus `H(F^n(T))=S^1` for every `n`.
For any angle `z`, the nested nonempty compact sets
`F^n(T) intersect H^{-1}(z)` show that `H(K)=S^1`.

The circle map has entropy `log(3)`: its `3^n` equally spaced preimages
of zero are separated by `2 pi/3` at some time before step `n`, and its
Lipschitz constant 3 gives the matching upper bound. Entropy cannot increase
under a continuous onto factor: pull back a finite open cover of the factor;
surjectivity preserves the minimum covering number of every iterated join.
This gives the stated lower bound on `K`.

Finally, `F(-w)=-F(w)` and `q(-w)=q(w)`, so
`H(-w)=H(w)+pi`. The map `exp(2iH(w))` descends to `K/{w~-w}` and still
factors its dynamics onto degree-minus-three angle multiplication.
The entropy bound is therefore not an artifact of retaining both ICA signs.

At `epsilon=0`, every point in `T` maps in one step onto `C=D=0` and the
attracting set is exactly the expanding circle. For positive `epsilon`, this
argument proves a factor, not a conjugacy or a solenoid classification.

## Numerical observations

The pilot tested all 80 values `epsilon=k/200`, `k=1,...,80`, from the
normalized vector `(0.83,0.29,-0.13,0.11)`. After 5,000 warmup steps,
36 values had positive leading tangent growth in every one of five
20,000-step windows: `k=1,...,35` and `k=41`. All pilot outcomes are retained.
Four prescribed continuations used 10,000 warmup steps and five million
observed steps, with an analytic tangent Jacobian and three-vector QR.

| epsilon | Mean tangent spectrum per iteration | Status |
| --- | --- | --- |
| 0.04 | `(1.096236, -3.170598, -3.575291)` | Inside the analytically certified parameter interval |
| 0.10 | `(1.010586, -1.229833, -1.870958)` | Numerical observation outside that interval |
| 0.15 | `(0.594634, -0.670846, -1.416793)` | Numerical observation outside that interval |
| 0.205 | `(0.168501, -0.489381, -1.481652)` | Numerical observation outside that interval |

Every 500,000-step window had positive leading growth in all four runs.
For each of `epsilon=0.04` and `0.15`, all 128 fresh starts passed the
positive-growth gate over 200,000 observed steps. The first 64 starts lie
in the stated tube; the other 64 are sampled on `S^3`. These are finite
cohorts, not a proof of a global basin.

Independent explicit-cubic calculations at 70 and 100 decimal digits each
retained positive leading growth over 20,000 observed steps for both portrait
parameters. The actual scikit-learn 1.9.1 public `FastICA` loop, using the
rounded 720-row data with `whiten=False`, also retained positive growth
through all 50,000 measured steps after 2,000 warmup steps. Its callable
nonlinearity delegated unchanged to the installed cubic function. Neither
first deflation unit stopped before the observation budget at `tol=1e-15`.
The largest step discrepancy against the polynomial model was `3.02e-13`.

The exact theorem concerns real radicals and exact coefficients. Rounded
data and machine arithmetic are separate numerical controls. At 1,024 random
test states across eight parameters, direct empirical updates and their
tangent derivatives agreed with the polynomial implementation to below
`9.8e-12` and `6.3e-12`, respectively; the latter is a relative tangent error.
An exact-zero-parameter circle check used 16,384 angles.

The portraits use 500,000 consecutive states each after 10,000 warmup steps,
projected onto `(A,B,C)`. They omit `D`; apparent crossings need not be
state-space intersections. Both panels use equal coordinate units and a
common scale, with no sign folding, trajectory replication or interpolation.
They show nested narrow strands at `epsilon=0.04` and broader structure at
`0.15`. Visual distinction from classical expanding-circle constructions is
not established.

## Prior art and scope

[Hyvarinen's FastICA paper](https://www.princeton.edu/~sswang/awake_cerebellar_ms_refs/hyvarinen99_ieee_trans_neural_netw.pdf)
defines the established fixed-point algorithm. This study changes its finite
input data, not that algorithm. [Wei's analysis, Section 2.1](https://arxiv.org/abs/1408.6693)
assumes a linear mixture of mutually independent sources; that statistical
assumption is not asserted for this constructed data set. The result is not
a counterexample to convergence theorems under those hypotheses.

The quartic representation makes the update an unshifted symmetric tensor
power iteration. Its failure to converge is already established:
[Kolda and Mayo, Section 3.2 and Example 3.5](https://arxiv.org/abs/1007.1267)
reproduce a fourth-order counterexample of Kofidis and Regalia.
[The original 2002 abstract](https://doi.org/10.1137/S0895479801387413)
states conditional convergence under convexity or concavity; its full text
was not recovered in this review. We claim neither the first nonconvergent
tensor iteration nor a new general chaos mechanism.

The contribution recorded here is an explicit 720-row realization, an exact
parameter interval with strict trapping, and the factor proof. Searches for
FastICA with tripling, solenoids, chaotic attractors and periodic orbits did
not locate this exact construction in the returned material. This is a bounded
search result, not a historical uniqueness certificate.

Transitivity of `K`, equality of `K` with a sampled orbit closure, uniform
hyperbolicity, a Smale--Williams conjugacy, and priority of the finite-data
construction remain unproved. The entropy theorem does not settle those
questions or establish the requested visual novelty.

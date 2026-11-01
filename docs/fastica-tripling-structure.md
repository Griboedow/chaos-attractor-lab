# Local noninjectivity at the FastICA attracting set

For every `0 < epsilon <= 1/25`, the attracting set `K` in the
[structured FastICA construction](fastica-tripling.md) contains a point at
which the **full three-dimensional map on the unit sphere is not locally
injective**. This excludes a conjugacy on an open neighborhood of `K` to a
map that is locally injective throughout its corresponding invariant set,
including the usual Smale--Williams solenoid embeddings.

This is an ambient-neighborhood obstruction. It does not prove that the
restricted map `F|K` is noninjective, identify the topology of `K`, or rule out
a conjugacy defined only on invariant sets. The argument also does not locate
or classify a particular Whitney fold.

[Exact checks and reproduction](../results/fastica-local-injectivity/README.md) ·
[Construction and entropy theorem](fastica-tripling.md) ·
[Main-page portraits](../README.md#structured-data-fastica)

## Argument

The strict trapping set

```math
T=\{(A,B,C,D)\in S^3:\sqrt{C^2+D^2}\leq\tfrac25\sqrt{A^2+B^2}\}
```

is a solid torus, hence connected. Each `F^n(T)` is nonempty, compact and
connected; they form a decreasing sequence. Their intersection `K` is
therefore connected. Every periodic orbit contained in `T` belongs to `K`,
since each point has a predecessor of every finite length in `T`.

Two symmetric families supply periodic points in `T` with opposite signs
of the oriented tangent determinant. If `F` were locally injective at every
point of `K`, invariance of domain would make it a local homeomorphism on
an open neighborhood of `K`. Its local orientation sign would then be a
locally constant function into `{+1,-1}`. On the connected set `K` that sign
must be constant, contradicting the two periodic points.

Invariance of domain is [Hatcher, Theorem 2B.3, p. 172](https://pi.math.cornell.edu/~hatcher/AT/AT.pdf#page=181).
The orientation argument uses the locally consistent generators of local
homology described on pp. 233--235 of the same book: a local homeomorphism
pulls back the target orientation, and its comparison with the source
orientation is locally constant. At a regular differentiable point this
comparison is the sign of the tangent Jacobian.

The remaining work is algebraic: construct the two periodic points,
place them inside `T`, and establish their opposite signs uniformly over
the stated parameter interval.

## Tangent determinant

Let `h=grad(P)/4` be the cubic numerator in the construction, and put
`J=Dh`. Its homogeneity gives `J(w)w=3h(w)`. At a unit point with
`h(w) != 0`, choose oriented orthonormal frames `(w,V)` and `(F(w),U)` in
four-dimensional space. Differentiating the normalization yields

```math
DF=\frac{(I-FF^{\mathsf T})J}{\|h\|},\qquad
\det(U^{\mathsf T}DFV)=\frac{\det J}{3\|h\|^4}.
```

For the determinant identity, express `J` in these two frames. Its first
column is `(3||h||,0,0,0)`, and its lower-right block is
`||h|| U^T DF V`. Both frame determinants are +1. Consequently the
oriented three-dimensional determinant has exactly the sign of `det J`.

## An axis fixed point

Consider the unnormalized vector `v_a=(1,0,t,0)`. The first numerator
component and the condition for a parallel image are

```math
\begin{aligned}
a_e(t)&=1+18et+36e^2t^2+18e^3t^3,\\
Q_a(t)&=18e^3t^4+(36e^2+18e^4)t^3
 +(18e-54e^3)t^2+(1-36e^2)t-6e,
\end{aligned}
```

where `e=epsilon` and `Q_a=t h_A-h_C`. The other two components vanish.
For `0<e<=E=1/25`,

```math
\begin{aligned}
Q_a(4e)/e &=2(2880e^6+720e^4+72e^2-1)<0,\\
Q_a(6e)/e &=216e^2(126e^4+27e^2+2)>0.
\end{aligned}
```

The first inequality follows by substituting `E` in its three positive
terms. Hence a root `t_a` exists in `(4e,6e)`. It is unique in that interval:
`Q_a'(t) >= 1-36E^2-648E^4 > 0` there.
Since `a_e(t_a)>0`, homogeneity gives

```math
w_a=\frac{(1,0,t_a,0)}{\sqrt{1+t_a^2}},\qquad F(w_a)=w_a.
```

Its transverse-to-planar coordinate ratio is `t_a<6E=0.24<0.4`, so
`w_a` lies in `T` and therefore in `K`. Direct differentiation gives

```math
\det J(1,0,t,0)
=-2916e^4(1+6e^2t^2)^3(4+18et+15e^2t^2)<0
```

at `t=t_a`. Normalizing the input multiplies this degree-eight
determinant by a positive factor and does not change its sign.

## A diagonal two-cycle

Now take `v_d=(1,1,t,t)`. The numerator satisfies `h_A=h_B` and `h_C=h_D`.
The corresponding expressions are

```math
\begin{aligned}
b_e(t)&=2(-1+9et+9e^2t^2+9e^3t^3),\\
Q_d(t)&=18e^3t^4+(18e^2+126e^4)t^3
 +(18e-54e^3)t^2-(2+18e^2)t-6e.
\end{aligned}
```

At the endpoints,

```math
\begin{aligned}
Q_d(-4e)/e&=2(1+180e^2-1008e^4-1728e^6)>0,\\
Q_d(-2e)/e&=-2(1-54e^2+180e^4+360e^6)<0.
\end{aligned}
```

For the first sign use `1-1008E^4-1728E^6>0`; for the second use
`1-54E^2>0`. A root `t_d` lies in `(-4e,-2e)` and is unique there because
`Q_d'(t) <= -2+1296E^4+6048E^6 < 0`. Throughout this interval,
`b_e(t) <= -2+288E^4 < 0`. Therefore

```math
w_d=\frac{(1,1,t_d,t_d)}{\sqrt{2+2t_d^2}},\qquad
F(w_d)=-w_d,\quad F(-w_d)=w_d.
```

The ratio is `|t_d|<4E=0.16<0.4`; both points lie in `T` and hence `K`.
On this symmetry plane,

```math
\det J(1,1,t,t)
=46656e^4(1+6et+3e^2t^2)^3(2-6et+15e^2t^2)>0.
```

Indeed, the first factor in parentheses is at least
`1-24E^2=601/625>0`, and the last is positive for `t<0`.

At the displayed parameter `e=1/25`, supplementary 90-digit calculations give:

| Family | t | Oriented tangent determinant |
| --- | ---: | ---: |
| Axis fixed point | 0.217923939862465 | -0.005728984012890 |
| Diagonal two-cycle | -0.113747837507034 | +0.003958244779071 |

The proof uses the uniform polynomial signs above, not these rounded values.
The [checker](../research/fastica_local_injectivity.py) verifies both
determinant factorizations, the endpoint identities, ten positive rational
margins, and supplementary oriented-frame and centered-difference controls.
It also independently checks nine uniform polynomial sign conditions using
exact Bernstein coefficients on `e^2 in [0,1/625]` and the relevant interval
of `t/e`. The reconstructed polynomial is checked symbolically; its Bernstein
basis functions are nonnegative and sum to one, so the smallest positive
coefficient supplies a strict lower bound on the whole rectangle.

## Consequences and limits of the comparison

The opposite signs and connectedness prove the announced failure of local
injectivity at some point of `K`. They also imply a zero tangent determinant
somewhere in `K`, but that weaker conclusion alone would not establish
noninjectivity: a singular differentiable map can still be locally one-to-one.

The obstruction survives ICA sign identification. The antipodal quotient
`S^3 -> RP^3` is a local homeomorphism. A locally injective descended map
near the image of a point would lift to a locally injective map near that
point, contradicting the result on `S^3`.

The standard solenoid comparison is precise. [Chen and Lin, Section I](https://web.ma.utexas.edu/mp_arc/c/14/14-41.pdf)
use `q(s,z)=(s^m, delta z+s/2)` with a contraction small enough to make
`q` an embedding of the solid torus. [Jiang, Wang, Zheng and Zhou,
Definition 3.8](https://www.impan.pl/shop/en/publication/transaction/download/product/89068)
define Smale solenoids through a homeomorphism and nested solid tori.
Such ambient maps are locally injective. No ambient homeomorphism can
conjugate our map on an open neighborhood of the entire `K` to such a map,
since local injectivity is preserved under conjugacy.

An intrinsic conjugacy `K -> Lambda` need not extend to a neighborhood.
The nearby pairs witnessing noninjectivity here need not both lie in `K`.
Thus this result leaves open injectivity of `F|K`, intrinsic solenoid
equivalence, transitivity, and equality between `K` and the numerical
portrait. It also does not distinguish the construction from all known
noninvertible systems. The repository already contains a different
[ambient obstruction for RAAR](raar-critical-membership.md).

The contribution is the explicit, uniform algebraic obstruction for this
finite-data FastICA realization. The topological tools, circle expansion,
tensor-power nonconvergence and noninvertible chaotic dynamics are prior
mathematics; historical uniqueness remains unresolved.

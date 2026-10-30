# Exact volume contraction in the motor model

The normalized motor field has an exact identity connecting phase-volume change to mechanical speed. In a parameter region containing both the certified switching saddle and the eight-row hyperchaos candidate, it proves exponential contraction of a smoothly weighted volume. Combined with the [compact global-attractor theorem](motor-global-attractor.md), this makes the global attractor a set of zero ordinary phase-space volume. It does not prove chaos, a fractal dimension, attraction to a switching network or scientific novelty.

## Field, constraint and divergence

Fix an integer $N\geq1$, $\rho>0$, $r>0$ and real $b$. Use the [source normalization](motor-rows-source-audit.md):

```math
V_j=-\rho X_j+r(s_j-\bar s),\qquad
c_j'=1-c_j-V_js_j,\qquad
s_j'=b-s_j+V_jc_j,\qquad X_j'=V_j,
```

with $\bar s=N^{-1}\sum_js_j$ and $\sum_jX_j=0$. Choose the $2N$ population coefficients and the first $N-1$ displacements as independent coordinates. The final displacement is minus the sum of the others. Their constant linear relation to an orthonormal coordinate system on the constraint plane changes volume by a constant only.

The diagonal derivatives give

```math
\operatorname{div}f=-2N-\rho(N-1)
+r\left(1-\frac1N\right)\sum_jc_j.
```

The forcing parameter $b$ does not appear in this instantaneous trace. The trace need not be negative at every state.

## An exact weighted-volume identity

Put $\vartheta=1-1/N$. The zero-sum mechanical constraint implies $\sum_jV_j=0$, hence

```math
\sum_jV_js_j=\frac1r\sum_jV_j^2
+\frac\rho r\sum_jX_jV_j.
```

Define the smooth scalar function and constant

```math
\Phi(y)=\vartheta\left(r\sum_jc_j+
\frac\rho2\sum_jX_j^2\right),\qquad
a_N=-2N+(N-1)(r-\rho).
```

Substituting the field, the mixed mechanical terms cancel exactly:

```math
\Phi'=r\vartheta\left(N-\sum_jc_j\right)
-\vartheta\sum_jV_j^2,\qquad
\operatorname{div}f+\Phi'=a_N-\vartheta\sum_jV_j^2.
```

Thus the divergence relative to the positive smooth density $e^{\Phi(y)}$ is at most $a_N$. Liouville's determinant identity supplies, for every forward trajectory,

```math
\det D\phi_t(y)=
\exp\left[a_Nt-\vartheta\int_0^t\sum_jV_j(u)^2\,du
-\Phi(\phi_t y)+\Phi(y)\right].
```

This determinant is positive. The same formula applies to the full $3N-1$ state variation; it is not a determinant restricted to a synchronized cluster subspace.

## Consequences and their precise scope

If

```math
r-\rho<\frac{2N}{N-1}\quad(N>1),
```

then $a_N<0$. On any compact positively invariant set $D$, the ordinary determinant obeys

```math
\det D\phi_t(y)\leq
\exp(\max_D\Phi-\min_D\Phi)\,e^{a_Nt}
\quad(y\in D).
```

Every bounded set enters such a $D$ by the global absorbing estimate. For $N=1$, the constraint fixes $X_0=0$ and the two population coordinates simply relax; the identity reduces to $\det D\phi_t=e^{-2t}$.

Let $\mathcal A_f$ be the compact strictly invariant global attractor. Its weighted volume $m_\Phi(\mathcal A_f)$ is finite. The change-of-variables formula and strict invariance imply

```math
m_\Phi(\mathcal A_f)=m_\Phi(\phi_t\mathcal A_f)
\leq e^{a_Nt}m_\Phi(\mathcal A_f).
```

When $a_N<0$, this forces $m_\Phi(\mathcal A_f)=0$. Since the density is positive and bounded above and below on the compact attractor, its ordinary Lebesgue measure is also zero. This statement is about the finite phase space; an infinite-dimensional PDE space has no corresponding canonical Lebesgue volume. Its attractor is the exact first-harmonic lift already described in the [PDE theorem](motor-global-attractor.md).

For any bounded trajectory where the Lyapunov exponents and time average below exist, dividing the determinant identity by time gives

```math
\sum_{k=1}^{3N-1}\lambda_k
=a_N-\vartheta\left\langle\sum_jV_j^2\right\rangle.
```

Equivalently, at finite time the sum of logarithms of all singular values equals the displayed determinant expression, without assuming an asymptotic limit. For a periodic orbit the endpoint term vanishes exactly. Any nonconstant periodic orbit with $N>1$ has a strictly positive mean square speed, so its sum of Floquet exponents is strictly below $a_N$. A negative sum is compatible with several positive individual exponents and does not establish an attracting chaotic set.

## The current source parameters

For $N=8$ and the fixed source value $r=1+(19/20)/\sqrt{32001/500000}$, the recorded exact parameter points all satisfy the strict volume condition:

| Parameter point | Exact $\rho$ | Upper bound $a_8$ on the exponent sum |
| --- | --- | --- |
| Attracting four-pair certificate | $3.3652183995790423$ | approximately $-6.27050646$ |
| Six-cluster saddle | $817/250$ | approximately $-5.58997767$ |
| Eight-row activity-path chaos candidate | $1291/440$ | approximately $-3.25261403$ |

The [exact algebra and interval parameter record](../results/motor-saddle-proof/volume-identity-recheck.json) gives the actual enclosures. These bounds concern every bounded trajectory at the named parameter point, independently of its observed dynamical regime. The last row does not upgrade the numerical candidate to a theorem about chaos. An initial report's manually transcribed candidate bound was corrected against this exact calculation; its original report and first algebra record are retained.

## Reproduction

```bash
python research/check_motor_volume_identity.py --output outputs/motor-volume-recheck.json
```

The [checker](../research/check_motor_volume_identity.py) differentiates the complete reduced polynomial field exactly for $N=1,2,3,4,8$, and verifies the constraint, mechanical cancellation, scalar derivative and weighted-divergence formulas. The derivation above proves the identities for general $N$; finite checks support its implementation. It also encloses the source parameter constants with Arb. No new trajectories are integrated, and no absence-of-prior-art claim is made.

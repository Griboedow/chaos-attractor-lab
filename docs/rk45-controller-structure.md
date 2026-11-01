# Structure and bounds of the RK45 controller candidates

The physical equation is the stable linear spring

\[
\dot y=Ay,\qquad A=\begin{pmatrix}0&1\\-1&-2\zeta\end{pmatrix},
\qquad 0<\zeta<1.
\]

Every exact solution tends to zero. The studied state is the numerical
algorithm's `(q/atol, p/atol, proposed h)`, not a physical three-dimensional
flow. The integrator is the existing Dormand–Prince 5(4) pair implemented by
[SciPy RK45](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html).
Its signed solution/error polynomials are

\[
R(z)=1+z+z^2/2+z^3/6+z^4/24+z^5/120+z^6/600,
\qquad E(z)=z^5(5z^2-39z+97)/120000.
\]

For a trial step, `v=R(hA)y`, `e=E(hA)y`, and

\[
s_i=\tau+\rho\max(|y_i|,|v_i|),\qquad
w=\sqrt{\tfrac12\sum_i(e_i/s_i)^2}.
\]

Here `tau=atol>0` and `rho=rtol>=0`. Acceptance requires `w<1`. The accepted
step uses the factor `min(10,0.9*w^(-1/5))`, additionally capped at one after
any rejected trial. Rejection reduces the trial step by
`max(0.2,0.9*w^(-1/5))` and retries the unchanged solution. Zero error selects
the maximum factor. The autonomous transcription omits the finite endpoint
and floating-point clock operations. The separate actual-library runs retain
the clock and all unmodified SciPy step operations.

The [priority audit](rk45-priority.md) credits the established polynomial-map,
controller-stability and error-norm framework. These formulas define the
specific experiment; they are not a new integrator.

## Symmetries and comparisons

Absolute tolerance scales out exactly when `u=y/tau`. Changing tolerance and
scaling the initial solution accordingly does not create a new system.
The complete map commutes with `(q,p,h)->(-q,-p,h)`; reflected shapes are
equivalent descriptions. The two D335 rosettes are visited by one trajectory.
A band exchange or its second iterate must not be counted as a second
attractor discovery.

The diagonalizing coordinate
`z=q-i*(p+zeta*q)/sqrt(1-zeta^2)` obeys
`z'=(-zeta+i*sqrt(1-zeta^2))*z`. Transforming the **complete** error rule
gives a conjugate algorithm. Replacing that rule by
`|error|/[tau+rho*max(|z_old|,|z_new|)]` changes the algorithm. This latter
control has exact phase symmetry, a two-dimensional amplitude/step quotient,
and a neutral phase direction for a bounded nonzero orbit. A periodic
quotient need not give a closed full-state orbit. Norm dependence is established
prior work, not a new principle inferred from this experiment.

## A sufficient accepted-solution bound

Let

\[
P=\begin{pmatrix}1&\zeta\\\zeta&1\end{pmatrix},\quad
r=\|y\|_P,\quad \lambda=-\zeta+i\sqrt{1-\zeta^2},\quad
a=\sqrt{2(1+\zeta)},\quad b=\sqrt{2(1+\zeta)/(1-\zeta)}.
\]

Suppose `|R(h lambda)|<=1` for `0<h<=H`, and
`|R(h lambda)/E(h lambda)|<=K` for every `h>=H`. If `K*b*rho<1`, every
sequence of accepted solution values, whatever the step proposal policy,
satisfies

\[
\boxed{\quad \|y_n\|_P\le\max\left(\|y_0\|_P,
\frac{K a\tau}{1-Kb\rho}\right).\quad}
\]

**Proof.** Direct algebra gives `A^T P+PA=-2*zeta*P` and, for any real
polynomial `f`, `f(hA)^T P f(hA)=|f(h lambda)|^2 P`. Thus
`r_new=|R|r` and `||e||_P=|E|r`. Acceptance implies

\[
\|e\|_2\le\sqrt2\max_i s_i
\le\sqrt2\tau+\sqrt2\rho\max(\|y\|_2,\|v\|_2),
\]

and the eigenvalues `1-zeta`, `1+zeta` of `P` give
`|E|r<=a*tau+b*rho*max(r,r_new)`. A step with `h<=H` cannot increase `r`.
For a larger step that does increase it, multiplication by `|R/E|<=K`
yields `r_new<=K*a*tau+K*b*rho*r_new`. Rearrangement gives the bound,
and induction completes the argument.

This elementary estimate bounds solution coordinates only. It does not
certify a full-state trapping region, chaos or historical novelty.
The broader use of error control to establish long-time numerical stability
is established; see [Stuart and Humphries (1994), Section 9](https://math.mcgill.ca/humphries/research/papers/1036054.pdf).
The calculation here supplies particular RK45 certificates, not a new general theory.

## Exact certificates used here

The error polynomial has no nonzero zero on the oscillator ray because

\[
|5(h\lambda)^2-39h\lambda+97|^2=
25h^4+390\zeta h^3+(551+1940\zeta^2)h^2+7566\zeta h+9409>0.
\]

In particular `|E(h lambda)|>=h^7/24000`. For a broad certificate use
`H=1/2`. The polynomial `(1-|R(h lambda)|^2)/h`, of degrees `(11,6)`
in `(h,zeta)`, has a nonnegative tensor Bernstein representation on
`[0,1/2] x [0,1]`. The audit saves all 84 exact rational coefficients and
reconstructs the polynomial exactly; this proves contraction throughout the
rectangle, rather than on a sampled grid. The triangle inequality then gives

\[
K_H=24000\sum_{k=0}^6 a_k H^{k-7}=5\,064\,880,
\qquad R(z)=\sum_{k=0}^6 a_k z^k.
\]

That constant is deliberately conservative. For selected parameter pairs,
interpreted as exact binary64 rationals, the audit obtains stronger certificates:

1. An exact Sturm count proves that `( |R(h lambda)|^2-1 )/h` has no zero
   on `[0,H]`; its value at zero is negative.
2. Every coefficient of
   `K^2*|E((H+t)lambda)|^2-|R((H+t)lambda)|^2` is nonnegative and its
   constant is positive. This proves the ratio bound for the entire
   unbounded interval `t>=0`.
3. The gain condition is checked as an exact rational inequality
   `(K*b*rho)^2<1`.

| Case | H | Squared gain | Normalized bound constant |
| --- | ---: | ---: | ---: |
| D335 | 3.3 | 2.89616e-18 | 3.88610 |
| D448 | 3.3 | 1.11928e-8 | 3.94639 |
| D2006 | 3.3 | 0.109630 | 10.13785 |
| D1890 | 2.2 | 15.0301 | Condition fails; no certificate claimed |

The last column is `K*a/(1-K*b*rho)`, so the bound is the larger of
the initial P-norm and this constant times `atol`. Failure of the sufficient
condition is not evidence of an unbounded orbit.

The reproducible audit is
[`audit_rk45_structure.py`](../research/audit_rk45_structure.py), with exact
coefficients in the frozen evidence packet. It derives `R` and `E` directly
from the rational seven-stage tableau before checking the inequalities.

## Limits

At exact zero solution, the error vanishes and the unrestricted proposed step
grows by a factor of ten each time. Consequently the full ideal map has no
compact global attractor over its entire state domain. A finite `max_step`
would define a different map. Bounded selected trajectories, positive finite
rates, their beautiful shape and the solution bound above do not prove
an asymptotic strange attractor or inequivalence to a known invariant set.

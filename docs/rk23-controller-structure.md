# Structure of the RK23 controller experiment

This note concerns the numerical algorithm for a stable linear oscillator. It
does not claim chaos in the differential equation, a new integrator, or a new
attractor class. Numerical evidence and the priority review are separate.

## Equation and accepted-step map

Let `0 < zeta < 1` and

\[
\dot q=p,\qquad \dot p=-q-2\zeta p,\qquad
A=\begin{pmatrix}0&1\\-1&-2\zeta\end{pmatrix}.
\]

Every exact solution tends to zero. With the Bogacki–Shampine 3(2) pair in
local extrapolation mode, the solution and embedded-error polynomials are

\[
R(z)=1+z+\tfrac12z^2+\tfrac16z^3,\qquad
E(z)=\tfrac1{48}z^3(1+z).
\]

For a trial step `h`, the new solution is `v=R(hA)y` and its error estimate
is `e=E(hA)y`. The experiment uses SciPy's componentwise scale

\[
s_i=\tau+\rho\max(|y_i|,|v_i|),\qquad
w=\sqrt{\frac12\sum_{i=1}^{2}(e_i/s_i)^2},
\]

where `tau=atol>0` and `rho=rtol>=0`. A step is accepted when `w<1`.
The next proposed step uses `min(10, 0.9*w^(-1/3))`; following a rejection
this factor is additionally capped at one. Rejected trials reduce their step
by `max(0.2, 0.9*w^(-1/3))` and retry the same solution state. Zero estimated
error selects the maximum factor. The transcription omits the finite final
time, finite-precision clock subtraction, and minimum-clock-step check;
unmodified-library runs retain these operations.

The error polynomial is the signed estimate in the inspected SciPy release.
Its sign does not affect the error norm. The fractions above define the
mathematical method; independently replaying its stages at high precision
tests that method rather than repeating rounded binary64 stage coefficients.

The general polynomial-map and step-control-stability framework is prior
work: see [Hall and Higham (1988)](https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P1.pdf)
and [Ranocha and Giesselmann (2024), Section 3](https://arxiv.org/html/2307.12677v3#S3).
The implementation policy is from
[SciPy RK23](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK23.html).

## Exact symmetries and coordinate dependence

The normalized solution `u=y/tau` removes the absolute tolerance from the
mathematical recurrence. Changing `tau` while fixing normalized initial data
therefore gives the same normalized dynamics. Starting from the same physical
displacement at a different tolerance changes that initial condition.

The real-coordinate map has the reflection symmetry
`(q,p,h) -> (-q,-p,h)`. It need not have continuous rotation symmetry:
the componentwise error scales and Euclidean RMS norm are taken in the
original axes. In contrast, a single complex equation with its complex
absolute-value error scale has exact phase-rotation symmetry. Its amplitude
and step form a two-dimensional quotient; adding the phase contributes a
neutral direction. Those rotated shapes are not independent new mechanisms.

Diagonalizing the physical oscillator does not by itself conjugate its
adaptive algorithm to that complex-scalar algorithm. The error norm and its
componentwise relative scales must also be transformed. This distinction is
already part of the step-control literature; it is not a novelty claim.

## An explicit bound for accepted solution values

The following deliberately conservative estimate certifies bounded solution
coordinates under a stated relative-tolerance condition. It is an auxiliary
calculation for this experiment, not a priority claim.

Put

\[
P=\begin{pmatrix}1&\zeta\\\zeta&1\end{pmatrix},\quad
\|y\|_P=\sqrt{y^TPy},\quad
\lambda=-\zeta+i\sqrt{1-\zeta^2},\quad
c=128\sqrt2.
\]

If `rho < (1-zeta)/c`, then **any sequence of accepted RK23 steps**, regardless
of the step proposal policy, satisfies

\[
\|y_n\|_P\le
\max\!\left(\|y_0\|_P,
\frac{c\tau\sqrt{1-\zeta}}{1-\zeta-c\rho}\right).
\]

This bounds positions and velocities, not the step-size coordinate.

**Proof.** Direct multiplication gives
`A^T P+PA=-2*zeta*P` and, for every real-coefficient polynomial `f`,

\[
f(hA)^TPf(hA)=|f(h\lambda)|^2P.
\]

Thus `r=||y||_P` changes to `r_new=|R(h lambda)| r`, and the error has
P-norm `|E(h lambda)| r`.

For `0<h<=1`, expansion yields

\[
|R(h\lambda)|^2-1=
-2\zeta h\left(1-\zeta h+\tfrac23\zeta^2h^2\right)
+\left(\tfrac23\zeta^2-\tfrac1{12}\right)h^4
-\tfrac16\zeta h^5+\tfrac1{36}h^6.
\]

Since `1-a+(2/3)a^2 >= 5/8`, `zeta^2 h^4 <= zeta h`, and `h^6<=h^4`,
the right side is at most `-(7/12) zeta h-h^4/18<0`. These steps reduce `r`.

For `h>=1`, let `omega=sqrt(1-zeta^2)`. The triangle inequality and
`|1+h lambda|^2=(h-zeta)^2+omega^2` give

\[
|R(h\lambda)|\le\tfrac83h^3,\qquad
|E(h\lambda)|\ge\tfrac1{48}h^3\omega,\qquad
\frac{|R|}{|E|}\le K:=\frac{128}{\omega}.
\]

Acceptance implies

\[
\|e\|_2\le\sqrt2\max_i s_i
\le\sqrt2\tau+\sqrt2\rho\max(\|y\|_2,\|v\|_2).
\]

Using the extremal eigenvalues `1-zeta` and `1+zeta` of `P`, this gives

\[
|E|r\le a\tau+b\rho\max(r,r_{\rm new}),\quad
a=\sqrt{2(1+\zeta)},\quad b=\sqrt{\frac{2(1+\zeta)}{1-\zeta}}.
\]

If the step increases `r`, multiplying by `|R|/|E|<=K` shows
`r_new <= Ka tau+Kb rho r_new`. The hypothesis is exactly `Kb rho<1`.
Rearrangement gives the displayed constant. Every step either decreases
the radius or lands below that constant; induction finishes the proof.

## A sharper certificate for O526

The broad sufficient condition above is deliberately conservative. Case O526
has `zeta=0.1717106344755739`, `rho=0.03929059404001495` and fails that condition.
A separate exact-rational calculation supplies a sharper bound for this
particular pair. Both numbers are interpreted as their exact binary64 values.
Put `H=23/10`. The polynomial
`Q(h)=(|R(h lambda)|^2-1)/h` has no roots on `[0,H]`, by an exact Sturm count,
and `Q(0)<0`. Hence all `0<h<=H` reduce the P-radius.

For `h>=H`, the triangle estimate can instead use

\[
K_H=\frac{48(H^{-3}+H^{-2}+(2H)^{-1}+1/6)}
{\sqrt{1-2\zeta H+H^2}}.
\]

The numerator decreases with `h` and the denominator increases for `h>=H>zeta`.
The earlier proof therefore gives
`||y_n||_P <= max(||y_0||_P, K_H*a*tau/(1-K_H*b*rho))`, provided
`K_H*b*rho<1`. The audit verifies its squared form as an exact rational
inequality and records the polynomial coefficients, rational parameters and
Sturm count. This is a solution-coordinate certificate for O526, not an
independent proof of hyperchaos. O37 and O689 satisfy the broad condition.

The reproducible symbolic audit is
[`research/audit_rk23_structure.py`](../research/audit_rk23_structure.py);
its frozen result is in the evidence packet's `source/` directory.

## Limits of these certificates

Bounded solution values and finite positive growth estimates are not a
proof of a strange attractor. The bound does not establish an invariant
trapping region in the full `(q,p,h)` state, physical invariant measures,
structural stability, or inequivalence to known attractors.

There is no compact global attractor for the ideal unrestricted-step map on
the whole state domain: at the exact solution `q=p=0`, the error is zero and
the proposed step grows by a factor of ten on every iteration. Reported
bounded candidate trajectories therefore require their own finite-data
step-size checks. Imposing a finite `max_step` would change the algorithm.

Persistent numerical oscillation near a stable physical equilibrium is
already described by
[Higham, Humphries and Wain (2000)](https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P35.pdf).
The present calculation must not be used to relabel that phenomenon as new.

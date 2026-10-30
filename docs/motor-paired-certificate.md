# An attracting four-pair motor cycle

The motor model has a **rigorously certified periodic orbit with four distinct synchronized pairs**, and that orbit attracts nearby trajectories in the full 23-dimensional eight-row system. The certificate gives one exact neutral flow multiplier and twenty-two stable multipliers. It turns the reference periodic state in the [postcritical numerical study](motor-postcritical.md) into a local computer-assisted theorem.

This is an attracting periodic orbit, not a strange attractor. Its existence and stability at one parameter do not prove the nonlinear bifurcation that creates it, connect it rigorously to the parent four-plus-four branch or establish scientific priority.

## Exact parameter and statement

The equations and source path are unchanged. Their exact inputs are

```math
\rho=3.3652183995790423,\qquad
r=1+\frac{19/20}{\sqrt{32001/500000}},\qquad b=1.
```

The finite decimal rho means an exact rational number. This differs slightly from the binary floating-point value used to initialize the numerical search. The exact algebraic gain is enclosed by ball arithmetic; no rounded decimal gain substitutes for it.

A nonconstant periodic solution exists with least positive period enclosed by

```math
T=6.69044151030481292225255\;\mathbin{\pm}\;9.19\times10^{-25}.
```

For that period, its full constrained-system Floquet spectrum has:

| Directions | Certified property |
| --- | --- |
| Flow direction | Exactly one multiplier equal to one |
| Other ten directions preserving the four pairs | Modulus strictly below 0.98 |
| Twelve directions splitting members of the pairs | Modulus strictly below 0.9 |

Thus the orbit is locally orbitally asymptotically stable. These are bounds on eigenvalue moduli, not a bound on one-step amplification in an arbitrary Euclidean norm. They do not determine the size of its basin.

All six comparisons between the four initial row states are separated by certified intervals. Each row state consists of `(c,s,X)`. The model is invariant under row permutations, and uniqueness of the flow prevents two distinct row states from coalescing later: equality at one time would imply equality along the entire orbit. The solution therefore has exactly four distinct row states, each repeated twice, throughout its motion.

The shooting root is locally unique in the stored phase-fixed box for the imposed quarter-period symmetry. No global uniqueness of periodic motion is claimed. A separate validated zero-count argument below establishes that the displayed period is the least one.

## Exact restriction and quarter-period closure

Write the four cluster states as `(c_a,s_a,X_a)`, with weights one quarter, and eliminate `X_3=-X_0-X_1-X_2`. The eleven independent coordinates are

```math
y=(c_0,c_1,c_2,c_3,s_0,s_1,s_2,s_3,X_0,X_1,X_2).
```

Their polynomial field is

```math
V_a=-\rho X_a+r(s_a-\bar s),\quad
\bar s=\frac14\sum_{b=0}^3s_b,\quad
c_a'=1-c_a-V_as_a,\quad
s_a'=1-s_a+V_ac_a,\quad X_a'=V_a.
```

Let `P` permute the row states in the order `(2,3,1,0)`. The constrained-coordinate matrix is exact and satisfies `P^4=I`. The field commutes with this permutation. We fix `X_0(0)=0` and solve eleven equations in the other ten initial coordinates and a quarter duration `L`:

```math
F(y_0,L)=\phi_L(y_0)-P y_0=0.
```

A validated Picard/Taylor flow, its eleven-dimensional fundamental matrix and a Krawczyk inclusion enclose a unique solution in the declared box. The quarter duration is positive and the initial section speed satisfies `V_0 > 1.1353`, excluding an equilibrium. Equivariance and flow uniqueness then give `phi_(4L)(y_0)=y_0`, so `T=4L` is a positive period.

The [448-bit recomputation](../results/motor-paired-proof/recheck-01.json) has Krawczyk contraction below `5.307e-22` and maximum relative image radius below `7.483e-19`. The initial ordinary Newton iterations are only a way to propose a center; none of their residual tolerances serves as an existence proof.

## Why this is the least period

The [least-period certificate](../results/motor-paired-proof/least-period.json) first revalidates the shooting root, then integrates only the eleven states with validated Picard tubes. Exact dyadic comparisons independently checked that its revalidated shooting box lies inside the spectral certificate's box, identifying the same locally unique root. Quarter symmetry lets it read the full `X_0` signal from the four row positions over one quarter, in the order `0,2,1,3`.

Every tube with a strict position sign excludes a zero. Maximal remaining windows have a strict velocity sign throughout; certified opposite signs at their boundaries therefore establish exactly one zero, by continuity and monotonicity. This proves existence as well as uniqueness of each counted crossing. The zero at total time zero is included, while the symmetry-related final zero at time T is excluded. All other quarter-boundary positions are separated from zero.

The count on `[0,T)` is exactly two upward and two downward zeros. A smaller state period would repeat the upward count an integer number of times, so the only alternative is `T/2`. But exact quarter symmetry gives `X_0(T/2)=X_1(0)>0`, whereas `X_0(0)=0`. This excludes `T/2`, establishing the least period T without inferring it from an approximate drawing or sampled return distance.

## All twenty-three tangent directions

Let `A=D phi_L(y_0)`. Differentiating the exact equivariance gives the quotient return matrix

```math
M_q=(P^{-1}A)^4.
```

An internal zero-sum contrast inside pair a leaves the mean response unchanged at first order and satisfies

```math
J_{\perp,a}=
\begin{pmatrix}
-1&-V_a-rs_a&\rho s_a\\
V_a&-1+rc_a&-\rho c_a\\
0&r&-\rho
\end{pmatrix}.
```

If `B_a` is its quarter-period fundamental matrix, row zero encounters the quarter blocks in the order `0,2,1,3`. Its complete return matrix is therefore `B_3 B_1 B_2 B_0`. The other three rows use the corresponding cyclic products. All four are enclosed and tested separately, although exact symmetry also makes their spectra equal.

The full tangent space is a direct sum of the eleven-dimensional quotient and four three-dimensional contrasts: `11+4*3=23`. This is an algebraic decomposition with a constant invertible embedding, not a projection of a computed leading vector. The controls independently assemble the original 23-dimensional Jacobian and check the block identity.

The autonomous nonconstant orbit supplies an exact factor `z-1` in the quotient characteristic polynomial. Synthetic division encloses the remaining degree-ten polynomial, using that exact analytic factor; a numerical eigenvalue close to one is not substituted for it. The quotient polynomial and the four transverse cubic polynomials then pass the following sufficient Schur test.

## Why the spectral test is rigorous

To test a radius R, first replace p(z) by p(Rz), then normalize its leading coefficient to one. For a real monic polynomial of degree n with constant coefficient a, put

```math
p^\#(z)=z^n p(1/z),\qquad
q(z)=\frac{p(z)-a p^\#(z)}{z}.
```

The interval code verifies `|a|<1` strictly and repeats the operation on q after normalization. The degree drops by one. Sufficiency follows by induction from the identity

```math
(1-a^2)p(z)=zq(z)+a q^\#(z).
```

If q has all its roots inside the unit disk, it has no root on the boundary. On that boundary `|q^#|=|q|`, so `|a q^#|<|zq|`. Rouché's theorem gives p the same number of interior roots as zq, namely n. The final nonzero constant starts the induction. This is the elementary real Schur recursion; it does not require approximating the roots or resolving their multiplicities numerically.

Every coefficient operation, normalization and strict inequality is enclosed with Arb. Replacing the normalized leading coefficient by exact one is a mathematical identity; other coefficient uncertainty is retained. A leading interval containing zero or any failed strict margin rejects the test. The final recomputation has smallest Schur margins above `0.0005678` for the radius-0.98 quotient test and `0.011996` for the radius-0.9 transverse tests.

## Physical and density-PDE interpretation

Each row satisfies

```math
\frac{d}{d\tau}(c_a^2+s_a^2)
=-2(c_a^2+s_a^2)+2(c_a+s_a).
```

At a maximum over a periodic orbit, Cauchy-Schwarz therefore gives `c_a^2+s_a^2 <= 2`. This bounds the whole periodic solution algebraically, rather than only at sampled drawing points. The source activity at this parameter is approximately `delta=0.2451657191`; the published source normalization and positivity checks apply as described in the [source audit](motor-rows-source-audit.md).

Under the sinusoidal-rate assumptions in the [density-PDE lifting result](motor-pde-lifting.md), the periodic orbit lifts exactly to an attracting periodic density solution. The omitted density components add only the stable return factor `exp(-T)` and do not change the unstable or neutral count.

There is also a finite-dimensional replication corollary. For any `N=4m`, replicate each of the four row states m times. The quotient field stays the same; there are `N-4` internal three-dimensional contrasts, all with the already certified stable transverse spectra. The replicated orbit has one neutral and `3N-2` stable multipliers. This deduction uses exact equal-cluster weights; it is not a new simulation for every N.

## Reproduce and inspect

Install the proof dependencies independently of the simulation/GPU stack, then use new output paths:

```bash
python -m pip install -r requirements-proof.txt
python research/validate_motor_paired_cycle.py --bits 384 --steps 256 --order 40 --radius 1e-36 --output outputs/paired-proof.json
python research/check_motor_paired_proof.py --output outputs/paired-controls.json
python research/validate_motor_paired_period.py --input results/motor-paired-proof/recheck-01.json --output outputs/paired-least-period.json
```

The [first successful certificate](../results/motor-paired-proof/attempt-02.json) uses 384 bits, 256 subdivisions and Taylor order 40. The [recomputation](../results/motor-paired-proof/recheck-01.json) changes these to 448 bits, 288 subdivisions and order 42, narrows the shooting box and adds the distinct-row and sharper-radius checks. The runs took approximately 155 and 143 seconds, respectively; the second uses the first root only as an initial center and recomputes the inclusion and spectrum.

The matrix version of the Taylor recurrence computes the eleven core jets once and propagates fundamental-matrix coefficients by

```math
M_{k+1}=\frac{L}{k+1}\sum_{j=0}^k J_j M_{k-j}.
```

This is the same polynomial recurrence as the generic augmented scalar implementation, executed through Arb matrix operations. The same routine is evaluated on the validated Picard tube to enclose the Taylor remainder. The `--generic-jets` option retains the scalar implementation for comparison.

The [successful generic-recurrence recheck](../results/motor-paired-proof/generic-recheck-02.json) independently executes every augmented scalar Taylor recurrence at 384 bits, 160 subdivisions and order 38. It revalidates the root, the four distinct cluster states and both spectral-radius bounds. This checks two implementations of the Taylor coefficient propagation, while retaining the same Picard/Taylor enclosure method and arithmetic library.

The [controls](../results/motor-paired-proof/analytic-controls-final.json) check 5,848 independently assembled field, Jacobian, permutation and full-space decomposition components at eight rational states. They also test the Schur recursion on prescribed stable real/complex roots and on unstable/boundary controls. At an exact equilibrium, 2,041 computed matrix Taylor coefficients enclose independently calculated exact rational values `(LJ)^k/k!`; the generic scalar recurrence agrees within its arithmetic bounds. A separate control ensures that a coefficient interval containing zero is retained unless it is exactly zero.

The unsuccessful first integration remains in the [proof directory](../results/motor-paired-proof/): its truncation enclosure was too wide for the chosen shooting box. A generic scalar-recurrence recheck with a wider shooting box certified periodic existence and all transverse blocks but left the quotient spectral test inconclusive; its periodic-only status is preserved. Before the successful run, coefficient omission was also changed to the exact-zero predicate `is_zero()`. An early exact-rational control used a ball-converting containment method, which can widen the target and reject a valid narrower enclosure; the retained failed control was corrected to compare exact dyadic endpoints with exact rationals. No failed attempt is treated as scientific counterevidence.

A read-only peer review checked the quarter-period compositions, Schur/Rouché argument and the matrix Taylor recurrence, including its Picard-tube remainder evaluation. The proof is auditable numerical code, not a proof-assistant formalization or external peer review. The custom integrator, Python-FLINT/Arb and their runtime form its computational basis. Source hashes and executed source snapshots accompany every integration result.

The implementation is [validate_motor_paired_cycle.py](../research/validate_motor_paired_cycle.py), with [independent controls](../research/check_motor_paired_proof.py) and a [separate least-period verifier](../research/validate_motor_paired_period.py). The surrounding [numerical branch continuation](motor-postcritical.md) and [priority audit](motor-followup-priority.md) provide context; their broader bifurcation interpretations remain separate from this theorem.

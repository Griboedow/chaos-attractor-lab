# A unique fourfold crossing on a four-pair cycle

A validated four-pair motor cycle has a **unique local transverse multiplier crossing with exactly four critical directions and nonzero crossing speed**. This certifies the candidate near rho = 3.2748747166 in the [numerical postcritical continuation](motor-postcritical.md). It is a local theorem; a rigorous continuation connecting it to the [earlier attracting four-pair cycle](motor-paired-certificate.md) has not been supplied.

The result does not establish a nonlinear flip bifurcation or the birth of another periodic branch. The fourfold degeneracy prevents a direct application of the generic simple-multiplier flip theorem. No chaotic attractor, heteroclinic connection, network attraction or scientific-priority claim follows from this certificate.

## Exact statement

Use the original 23-dimensional relative-coordinate system for eight rows, the normalization in the [source audit](motor-rows-source-audit.md), and the fixed exact parameters

```math
r=1+\frac{19/20}{\sqrt{32001/500000}},\qquad b=1.
```

The following finite decimal endpoints mean exact rational numbers:

```math
I_4=[3.2748747165675179794363952840257987573543731832880148565086386329,
\quad 3.2748747165675179794363952840257987573543731832880148765086386329].
```

For every rho in this interval, the stored common shooting box contains a locally unique periodic solution with four distinct cluster states, each occupied by two rows, satisfying the quarter-period permutation closure. Its phase section is transverse. Exactly one rho in the interior has a transverse multiplier equal to minus one. At that parameter:

- Minus one has algebraic and geometric multiplicity exactly four in the full 23-dimensional system.
- Eighteen multipliers lie strictly inside the unit circle.
- The remaining multiplier is the exact, simple flow multiplier one.
- The critical multiplier has the outward crossing-speed bounds

```math
18.827178 < \frac{d\mu_*}{d\rho} < 18.827180.
```

The lower endpoint has four unstable, one neutral and eighteen stable multipliers. The upper endpoint has no unstable, one neutral and twenty-two stable multipliers. Increasing source activity decreases rho, so increasing activity crosses this boundary toward instability.

Throughout the common family, the three simple real roots of each transverse cubic lie in `(-1.001,-0.999)`, `(-0.0021,-0.0018)` and `(0.320,0.323)`. All ten nonflow quotient multipliers have modulus strictly below 0.95. These separations rule out additional critical directions.

The least positive period is approximately `T = 6.81054204975892796672`, and the corresponding source activity is approximately `delta = 0.31481389327985162394`, using `delta=(r-1)^2/rho^2-1`. Exact dyadic enclosures in the evidence files are authoritative. The many parameter digits describe the exact mathematical interpretation of the source constants, not experimental precision. The very narrow family is chosen to control interval wrapping.

## Closure and exact multiplicity

The quotient coordinates are the four c values, four s values and three independent X values, with `X3=-X0-X1-X2`. Set `X0(0)=0` and solve eleven equations for ten free initial coordinates and the quarter time L:

```math
F(q,\rho)=\phi_L(y_0)-P y_0=0,
\qquad P=(2,3,1,0),\qquad P^4=I,\qquad T=4L.
```

The model is exactly equivariant under this equal-weight row permutation. Therefore the closure supplies a complete periodic orbit. Pairwise c/s/X differences certify four distinct row states; uniqueness of the row equations prevents two states from coalescing later. The positive section velocity makes the orbit nonconstant.

The full tangent system splits exactly into the eleven-dimensional quotient and four independent three-dimensional within-pair contrasts. With `A=D phi_L(y0)`, the quotient return is `(P^{-1}A)^4`. If Bj is the quarter fundamental matrix for the internal contrast of row j, the row-zero return is

```math
M_0=B_3 B_1 B_2 B_0.
```

The other row paths are cyclic rotations of `0,2,1,3`. Every fundamental matrix is invertible. In particular,

```math
M_2=B_0 M_0 B_0^{-1}.
```

The analogous cyclic identities make all four characteristic polynomials, and hence all four functions `g(rho)=det(I+Mj(rho))`, exactly equal. This identity is essential: separate interval brackets alone would not prove that four critical roots occur at the same parameter. Three disjoint sign-changing intervals account for the entire degree-three polynomial and make each of its roots simple. The exact tangent-space decomposition then gives precisely four independent copies of the critical root.

## Analytic parameter sensitivity

A uniform Krawczyk contraction establishes the common shooting family. Its center residual uses a parameter mean-value enclosure, while the complete interval shooting Jacobian remains enclosed. All integrations use validated Picard tubes and interval Taylor remainders. A rigorous ball-matrix solve gives the implicit derivative

```math
q'(\rho)=-[D_qF]^{-1}\partial_\rho F.
```

For normalized time u, let `x_u=L H(x,rho)` and put `alpha=L'/L`. The total state sensitivity w satisfies

```math
w_u=L\bigl(D_xH\,w+\partial_\rho H+\alpha H\bigr).
```

The initial w is the corresponding derivative of the shooting coordinates. The fixed-shooting partial derivative is integrated separately with zero initial w and alpha zero. No finite difference of integrated trajectories supplies a derivative bound.

For each transverse block, a six-dimensional triangular variational system propagates its fundamental matrix B and total derivative D simultaneously:

```math
\begin{pmatrix}B&0\\D&B\end{pmatrix}_u
=L\begin{pmatrix}J&0\\J'+\alpha J&J\end{pmatrix}
\begin{pmatrix}B&0\\D&B\end{pmatrix}.
```

Here J' includes the analytic parameter derivative and the state sensitivity. The initial triangular matrix is the identity. The full return derivative is then assembled with the product rule in the certified quarter order. A determinant derivative computed by multilinearity remains valid when `I+M` is singular; no inverse is used.

For `p(lambda,rho)=det(lambda I-M(rho))`, the resulting identities are

```math
g(\rho)=-p(-1,\rho),\qquad
\mu_*'(\rho_*)=\frac{g'(\rho_*)}{\partial_\lambda p(-1,\rho_*)}.
```

All four block enclosures give `g'` near `24.8310533909907731701954`, with displayed radii below `4.2e-23`, and strictly positive critical speed. Fresh endpoint shooting boxes lie strictly inside the common family box. Their determinant signs are strictly negative and positive, respectively, near `-2.4831e-52` and `+2.4831e-52`. Continuity gives a crossing, strict positivity of g' makes it unique, and the root separation gives its exact multiplicity.

The quotient uses its exact autonomous factor `lambda-1`, followed by the interval Schur recursion described in the [paired-cycle certificate](motor-paired-certificate.md). A numerical eigenvalue near one is not substituted for this analytic factor.

## Least period on the whole family

The separate [least-period certificate](../results/motor-paired-crossing-proof/least-period.json) freshly validates a larger common shooting box and checks exact dyadic containment of the crossing family's box and parameter interval. Local uniqueness therefore identifies the same periodic solutions. State-only validated Picard tubes cover the complete quarter orbit at every parameter in the common interval.

Under the quarter permutation, total-time X0 follows quotient rows `0,2,1,3`. All possible zero windows are resolved by a strict derivative sign together with certified endpoint signs. Opposite signs prove an interior zero exists and strict monotonicity makes it unique. The initial upward zero at time zero is counted; the corresponding final zero at time T is excluded. This gives exactly two upward and two downward zeros in `[0,T)`.

If a nonconstant orbit had a smaller positive state period, its minimal period would divide T by an integer, and the upward-zero count would repeat that many times. Only T/2 could therefore remain. Exact quarter symmetry gives `X0(T/2)=X1(0)>0`, whereas `X0(0)=0`, excluding T/2. Consequently T is the least period throughout the certified crossing family. This does not supply the nonlinear hypotheses for a flip bifurcation.

## Evidence and reproduction

The [successful certificate](../results/motor-paired-crossing-proof/attempt-02.json) uses 512-bit Arb arithmetic, 288 subdivisions and Taylor order 64. Its common-box Krawczyk contraction is below `1.870e-35`, and the maximum inclusion ratio is below `0.002123`. Every strict gate passes. The largest recorded Picard contraction is below 0.560. The run took approximately 2,691 seconds, including three center proposals; the proposals are guesses and are not counted as proof.

The [changed-precision recomputation](../results/motor-paired-crossing-proof/recheck-01.json) uses 576 bits, 320 subdivisions and Taylor order 68. It takes the earlier center only as a guess and freshly recomputes the uniform inclusion, analytic sensitivities, complete spectrum, endpoint inclusions and determinant signs. All gates pass again; the run took approximately 2,179 seconds. The [exact link audit](../results/motor-paired-crossing-proof/proof-links.json) verifies that both runs use the same exact finite parameter interval and that both shooting boxes and their endpoint parameter enclosures fit inside the separate least-period anchor. It also verifies the input hashes, identical shared sources and all seventeen archived integration-source snapshots. Thus the least-period conclusion applies to both complete crossing certificates.

The [independent controls](../results/motor-paired-crossing-proof/analytic-controls-final.json) pass:

- 1,280 independently assembled state, parameter, duration and triangular-field components;
- 4,018 nonlinear jet coefficients compared with the full generic sparse recurrence;
- 3,445 exact rational constant-block Taylor coefficients;
- 101 exact rational product-rule, determinant-derivative and cyclic-invariance checks, including a singular `I+M` example with known nonzero derivative.

The [first attempt](../results/motor-paired-crossing-proof/attempt-01.json) remains archived. Its sufficient Picard contraction bound exceeded one during a center proposal. The correction uses separate positive diagonal norm weights for state sensitivities, fundamental matrices and derivative blocks. All coordinates and interval uncertainties remain present; only the norm for the contraction test changes. This failed attempt is not evidence against the orbit or crossing.

With the proof dependencies installed, choose unused output paths:

```bash
python research/validate_motor_paired_crossing.py --output outputs/fourfold-crossing.json
python research/check_motor_paired_crossing.py --output outputs/fourfold-controls.json
python research/validate_motor_paired_crossing.py --center-input results/motor-paired-crossing-proof/attempt-02.json --refine-iterations 0 --bits 576 --steps 320 --order 68 --output outputs/fourfold-recheck.json
python research/validate_motor_paired_crossing_period.py --input results/motor-paired-crossing-proof/attempt-02.json --output outputs/fourfold-least-period.json
python research/check_motor_paired_crossing_links.py --first results/motor-paired-crossing-proof/attempt-02.json --recheck results/motor-paired-crossing-proof/recheck-01.json --period results/motor-paired-crossing-proof/least-period.json --output outputs/fourfold-proof-links.json
```

Every integration result retains the executed source snapshots and SHA-256 hashes. The implementation is [validate_motor_paired_crossing.py](../research/validate_motor_paired_crossing.py), the controls are [check_motor_paired_crossing.py](../research/check_motor_paired_crossing.py), and the separate period verifier is [validate_motor_paired_crossing_period.py](../research/validate_motor_paired_crossing_period.py). Read-only review checked the normalized sensitivity equations, triangular blocks, product order, determinant sign relation and uniform parameter inclusion. The period proof reuses the previously reviewed zero-window argument and classification. This is auditable numerical proof code, not proof-assistant formalization or external peer review; the custom interval integrator and Python-FLINT/Arb remain its computational basis.

The [packet manifest](../results/motor-paired-crossing-proof/manifest.json) inventories the report, live source files and complete evidence subtree with raw-byte SHA-256 hashes. The manifest itself is excluded from its own inventory. The containment and provenance checker is [check_motor_paired_crossing_links.py](../research/check_motor_paired_crossing_links.py).

## Replication and density lift

For any `N=4m` with m at least two, replicate each of the four row states m times. The quotient field is unchanged, and the `N-4` internal contrasts have the same certified transverse spectrum. At the crossing there are exactly `N-4` copies of minus one, `2N+2` stable multipliers and one flow multiplier. The lower endpoint has `N-4` unstable directions; the upper endpoint has `3N-2` stable directions and one neutral direction. For N=4 there are no internal contrasts, so this transverse crossing is absent from that quotient system.

Under the sinusoidal-rate assumptions of the [density-PDE lifting theorem](motor-pde-lifting.md), omitted density components contribute only the stable return factor `exp(-T)`. They do not change the critical or unstable counts. This local spectral conclusion does not certify the intervening global branch, any connection to the six-cluster saddle, or an attracting heteroclinic network.

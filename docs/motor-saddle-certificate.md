# A certified six-cluster periodic saddle

The published eight-row motor model has a **rigorously certified periodic saddle with one expanding direction** at the exact parameter `rho=3.268`. Four rows are individual clusters and the remaining four form two synchronized pairs. The full constrained system has 23 state coordinates; the certificate accounts for all of them.

This establishes the periodic states underlying the [numerical switching study](motor-postcritical.md). It does not establish the infinite-time connections between those states, attraction of a network or a chaotic attractor.

## Exact parameters and spectral statement

The unchanged polynomial motor equations use

```math
\rho=\frac{817}{250},\qquad
r=1+\frac{19/20}{\sqrt{32001/500000}},\qquad b=1.
```

The source activity increment is exactly

```math
\delta=\frac{451250}{32001\rho^2}-1
       =\frac{18955151}{59169849}
       \simeq0.3203515188960513.
```

With zero-based row labels, the clusters are `(0),(1),(2,3),(4),(5),(6,7)`. A positive period is enclosed by

```math
T=13.6331457678306460232747\;\mathbin{\pm}\;9.71\times10^{-24}.
```

The [combined certificate](../results/motor-saddle-proof/similarity-certificate.json) proves:

| Tangent directions | Certified property |
| --- | --- |
| Flow direction in the 17D cluster quotient | One multiplier exactly equal to one |
| Other sixteen quotient directions | Multiplier moduli below 0.8 |
| Three contrasts splitting pair `(2,3)` | Multiplier moduli below 0.9 |
| Contrasts splitting pair `(6,7)` | One simple real multiplier in `(1.12,1.13)`; two moduli below 0.9 |

Thus the complete count is **one unstable, one neutral and twenty-one stable**. The stable bounds concern eigenvalue moduli, not the operator norm in arbitrary coordinates. A periodic saddle is not an attracting periodic orbit in the full space.

All fifteen comparisons between the six initial `(c,s,X)` cluster triples are separated by intervals. Uniqueness of the flow prevents two distinct row states from later becoming identical, so this cluster partition holds throughout the orbit. The initial section velocity satisfies `X_0' > 0.4027`, excluding an equilibrium. The result is locally unique in the stored phase-fixed shooting box with its imposed symmetry; no global uniqueness is claimed.

## Weighted restriction and multiple shooting

The six cluster weights are

```math
w=\frac18(1,1,2,1,1,2),\qquad
\bar s=\sum_{a=0}^{5}w_a s_a,\qquad
X_5=-\frac{\sum_{a=0}^{4}w_aX_a}{w_5}.
```

The quotient has six `c`, six `s` and five independent `X` coordinates. It obeys `V_a=-rho*X_a+r*(s_a-mean_s)`, `c_a'=1-c_a-V_a*s_a`, `s_a'=1-s_a+V_a*c_a`, and `X_a'=V_a`. The unequal weights are essential.

Let `P` permute the cluster states in the order `(1,0,2,4,3,5)`. It preserves the weights, commutes with the field and satisfies `P^2=I`. The orbit is closed by

```math
\phi_{T/2}(y_0)=P y_0.
```

A single long interval integration was too wide. The successful calculation instead divides this half period into four equal durations `L=T/8`. It solves

```math
\phi_L(y_i)-y_{i+1}=0\quad(i=0,1,2),\qquad
\phi_L(y_3)-P y_0=0.
```

Fixing `X_0(0)=0` leaves `16+3*17+1=68` unknowns and equations. Every segment has its own validated Picard tube and Taylor remainder. The Krawczyk test uses the full 68-dimensional derivative, including the common duration column. Exact matching at the root joins the four segments into a single trajectory; equivariance supplies the second half.

The [first successful existence run](../results/motor-saddle-proof/attempt-02.json) uses 448-bit arithmetic, 128 subdivisions **per segment**, Taylor order 44 and coordinate box radius `1e-35`. Its maximum relative Krawczyk image radius is below 0.013636, and its contraction is below `1.548e-17`. The ordinary Newton residuals only select the center and do not prove existence.

## Monodromy and stable-root certification

If `M_i` is the quotient fundamental matrix of segment i, then

```math
H=M_3M_2M_1M_0,\qquad M_q=(PH)^2.
```

Each internal pair contrast has the three-dimensional Jacobian derived in the [cluster decomposition](motor-followup-priority.md). Both paired clusters are fixed by `P`, so each contrast return matrix is the square of its own four-segment product. The exact decomposition is `17+3+3=23`.

The expanding cubic has opposite strict signs at 1.12 and 1.13, enclosing a real root. Synthetic division by an interval containing that exact root gives a quadratic which passes a strict radius-0.9 Schur test. Its two roots therefore lie inside the disk, and the expanding root is simple. The other contrast cubic passes the same disk test directly.

The quotient required an additional algebraic step. The original characteristic-polynomial enclosure was too broad for the Schur test, even though the matrix enclosure was valid. That inconclusive result remains in `attempt-02.json`; it was not changed into a success.

An ordinary real-Schur computation proposes a matrix `Q`. Its entries are then stored as **exact dyadic rational numbers** in the [basis record](../results/motor-saddle-proof/spectral-basis.json). Arb verifies that its determinant excludes zero and encloses its inverse. The proof evaluates

```math
\widetilde M=Q^{-1}M_qQ
```

with every original matrix radius retained. No approximate eigenvalue is accepted, and `Q^T` is not substituted for the inverse. Any exactly specified nonsingular basis would be valid; the numerical proposal is useful only because it yields tighter bounds.

Exact similarity preserves the characteristic polynomial. The certified nonconstant autonomous orbit supplies a factor `z-1`; interval synthetic division encloses the remaining degree-sixteen polynomial. Its radius-0.8 Schur recursion has every strict margin above 0.189. The [Schur argument](motor-paired-certificate.md) uses Rouché's theorem and does not require numerical root separation. The stored flow-factor remainder contains zero, but that check alone is not the reason the exact factor exists.

The combined certificate reuses the validated matrices from `attempt-02.json` and clearly marks the similarity step as postprocessing. It does not describe that step as a repeated integration. The [postprocessor](../research/certify_motor_saddle_spectrum.py) and the optional `--spectral-basis` integration path both preserve the interval uncertainty.

A [fresh interval integration](../results/motor-saddle-proof/recheck-01.json) also passes every existence and spectral gate at 512 bits, 160 subdivisions per segment and Taylor order 46. It starts from the archived center, with no new floating-point Newton iteration. A [second algebraic postprocessing record](../results/motor-saddle-proof/recheck-postprocess-01.json) preserves both the older integrator-time helper and the subsequently executed helper under different snapshot names; every input uncertainty and source hash remains traceable.

## Least periods

The [least-period certificate](../results/motor-saddle-proof/least-period.json) revalidates a shooting box contained in the spectral certificate's box, linking its trajectory to the same locally unique root. Validated state tubes show exactly four upward and four downward zeros of `X_0` on `[0,T)`. The initial upward zero is counted once and the final repeated zero is excluded. Any smaller state period must therefore divide this count, leaving only `T/2` or `T/4` as candidates. But `X_0(T/2)=X_1(0)` lies near `-0.1950599752`, whereas `X_0(0)=0`. This excludes both possibilities and proves that **T is the least period of the complete state**.

The synchronized row `(2,3)` is fixed by the half-period permutation. Its displacement has two upward and two downward zeros during that half period, while `X_2(T/4)-X_2(0)` lies near `0.7891021643` and excludes zero. Hence its displacement and complete row state have least period **T/2**. These statements use monotonic zero windows enclosed by Picard tubes, not sampled sign changes.

## Leading expansion and contraction rates

The [rate certificate](../results/motor-saddle-proof/leading-rates.json) refines the two transverse cubic polynomials from the fresh interval integration. Strict rational endpoint signs enclose a positive real root in each block; division by its enclosing interval and a radius-0.11 Schur test locate both remaining roots strictly inside that disk. The stable root exceeds 0.8, so it also dominates every nonflow quotient multiplier. The certified leading rates per normalized time are

```math
e=\frac{\log\mu_u}{T}
 =0.00869626779094932\;\mathbin{\pm}\;3.79\times10^{-18},
\qquad
c=-\frac{\log\mu_s}{T}
 =0.0163144859413539\;\mathbin{\pm}\;3.58\times10^{-17},
```

and

```math
\frac ce=1.87603306769524\;\mathbin{\pm}\;6.41\times10^{-15}>1.
```

These are local Floquet rates. Their ratio helps interpret the growing numerical residence times, but does not by itself prove a connecting orbit, attraction of a network or an exact recurrence for switching times. The separately published forecast used earlier floating-point estimates; that historical prediction is not replaced by these later bounds.

A [postprocessor recheck](../results/motor-saddle-proof/leading-rates-recheck.json) gives tighter overlapping intervals. Five [known-cubic controls](../results/motor-saddle-proof/leading-rate-controls-02.json) pass, including real and complex stable remainders and rejection of a second large or repeated large root. The first control attempt stopped with inconclusive coarse intervals when an exact test root equalled the bisection midpoint. The revised search probes quarter endpoints in that case; the original failed control and already successful motor-rate certificate remain unchanged, and no acceptance condition was relaxed.

## Consequences and limits

Within the invariant subspace `q_6=q_7`, the expanding contrast is removed. The orbit is then locally orbitally attracting in twenty dimensions, with nineteen stable directions and the neutral flow direction. In the different subspace `q_2=q_3`, it remains a saddle. A row permutation exchanges the roles needed for the numerical source-to-destination experiment. The [connection strategy](motor-network-proof-strategy.md) explains why these local facts still leave a global existence problem.

The exact half-period symmetry and the six distinct clusters supply the hypotheses of the [finite-group lemma](motor-network-combinatorics.md). Consequently, the generated group has **eight distinct periodic-orbit images**. Their candidate graph has sixteen edges from the two observed numerical departures and their symmetry images. The distinction is essential: the node orbits are certified; the heteroclinic edges remain unproved. These are symmetry copies of one saddle, not eight newly discovered attractors.

The [full isotropy extension](motor-saddle-isotropy.md) combines the stable-versus-unstable pair spectra with the row's least period. It proves that the complete orbital stabilizer in `S_8` has order eight, giving **5,040 distinct labeled images** of this same saddle. Extending only the two specified departure rules gives 630 formal eight-node components; this does not exclude other dynamical connections.

Under the [density-PDE lifting assumptions](motor-pde-lifting.md), the same orbit lifts to the sinusoidal density model. Higher density modes contribute only stable return factors `exp(-T)`. The full PDE therefore retains the one unstable direction and the flow direction, with an infinite stable complement. This statement does not cover arbitrary transition rates or a bending flagellum coupled to fluid.

## Reproduce and audit

Use the optional proof environment and new output paths:

```bash
python -m pip install -r requirements-proof.txt
python research/check_motor_saddle_proof.py --output outputs/saddle-controls.json
python research/validate_motor_saddle_cycle.py --spectral-basis results/motor-saddle-proof/spectral-basis.json --output outputs/saddle-proof.json
python research/certify_motor_saddle_spectrum.py --input results/motor-saddle-proof/attempt-02.json --basis results/motor-saddle-proof/spectral-basis.json --output outputs/saddle-postprocess.json
python research/validate_motor_saddle_period.py --input results/motor-saddle-proof/similarity-certificate.json --output outputs/saddle-period.json
python research/certify_motor_saddle_rates.py --input results/motor-saddle-proof/recheck-01.json --output outputs/saddle-rates.json
```

The postprocessor checks an archived matrix enclosure; the integration command constructs a new existence and spectrum certificate. The period command revalidates existence and counts the zero windows. None requires the simulation/GPU environment. The approximate seed states and proposed similarity basis are already included and do not supply unverified accuracy assumptions.

[Independent algebra controls](../results/motor-saddle-proof/analytic-controls-01.json) check 7,144 field, Jacobian, weighted restriction, permutation and full-space block components at eight rational states. Another 3,991 matrix-jet coefficients enclose exact rational equilibrium references. These checks support the implementation; they are not identities proved by a finite sample alone.

A [separate review control](../results/motor-saddle-proof/independent-review-controls.json) compares nonconstant scalar and matrix Taylor recurrences, tests all 68 multiple-shooting residuals and 4,624 Jacobian entries on an exactly soluble affine field, and distinguishes correct from reversed products using noncommuting matrices. Source snapshots and hashes accompany successful and unsuccessful attempts.

Eight [similarity controls](../results/motor-saddle-proof/independent-similarity-controls/analysis.json) check a nonorthogonal exact basis, stable real and complex spectra, retained interval uncertainty, and rejection of an unstable root, an extra unit root, excessive uncertainty, singularity and nondyadic entries. Seven [zero-window controls](../results/motor-saddle-proof/zero-window-controls.json) exercise boundary ownership and rejection cases. These are implementation checks in addition to the orbit certificates.

The first single-shot attempt failed a Picard tube test. The first multiple-shooting run certified existence and both contrast blocks but left quotient stability inconclusive before the exact similarity step. Both histories are retained. This is a computer-assisted proof built on custom interval code and Python-FLINT/Arb, not a proof-assistant formalization or external peer review. The [FLINT documentation](https://flintlib.org/doc/arb_mat.html) describes the enclosing matrix operations used here; approximate matrix operations are not used in the certificate.

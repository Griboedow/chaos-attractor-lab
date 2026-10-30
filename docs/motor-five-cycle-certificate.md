# An attracting five-cluster motor cycle

An interval calculation certifies an **attracting periodic orbit with five distinct cluster states** in the full eight-row motor model. Its cluster multiplicities are `2+2+1+2+1`. Exactly one Floquet multiplier is the neutral flow multiplier; all twenty-two others have modulus strictly below 0.13.

This is a periodic attractor, not a chaotic attractor. The calculation does not establish historical novelty, a global branch connection, or the fraction of initial conditions in its basin. It certifies a positive return period; its least period has not been proved.

## Exact statement

Use the 23-dimensional relative-coordinate motor system from the [source audit](motor-rows-source-audit.md), with exact parameters

```math
\rho=\frac{129}{55},\qquad
r=1+\frac{19/20}{\sqrt{32001/500000}},\qquad b=1.
```

The stored multiple-shooting box contains a nonconstant periodic solution with a positive period enclosed by the outward readable bounds

```math
25.0016843758970971 < T < 25.0016843758970972.
```

Its five cluster states are pairwise distinct. One labeling, in zero-based row indices, is

```text
{0,3}, {1,2}, {4}, {5,6}, {7}.
```

The full tangent space has one simple multiplier one and twenty-two multipliers satisfying `|mu|<0.13`. Consequently the periodic orbit is locally orbitally exponentially attracting. This is a local stability theorem; no validated trajectory from a particular distant initial condition or quantitative basin boundary is supplied.

The numerical seed came from case 16 of the separately archived [four-probe destination study](../results/motor-activity-continuation/probe-destinations/README.md). That study followed all four physical perturbations of one earlier saddle endpoint. Three returned to symmetry-related copies of the old saddle; the fourth approached this different periodic state. Full DOP853 and Radau variational matrices agreed before interval validation was attempted. Those observations provide context and initial guesses, not premises of the interval theorem. The exact rational rho above also differs slightly from its binary floating-point grid seed; the certificate recomputes the root at the stated exact parameter.

## Weighted symmetry and multiple shooting

Write one c, s and X value for each cluster, with weights

```math
w=(2,2,1,2,1)/8,\qquad
\sum_{j=0}^{4}w_jX_j=0.
```

Eliminating X4 gives fourteen quotient coordinates. The weighted field is

```math
V_j=-\rho X_j+r\left(s_j-\sum_k w_k s_k\right),\qquad
c_j'=1-c_j-V_js_j,\quad s_j'=1-s_j+V_jc_j,\quad X_j'=V_j.
```

The permutation `P=(3,0,4,1,2)` preserves the weights and has order six. It cycles the three paired clusters and swaps the two singleton clusters. The validator solves the relative closure `phi_L(y0)=P y0` with `L=3 ell`, using three separate shooting segments:

```math
\phi_\ell(y_0)-y_1=0,\qquad
\phi_\ell(y_1)-y_2=0,\qquad
\phi_\ell(y_2)-P y_0=0.
```

The phase condition is `X0(0)=0`. Thirteen free coordinates in y0, fourteen in each of y1 and y2, and ell give forty-two unknowns and equations. A strict Krawczyk inclusion validates the root. The segment duration and initial section velocity are strictly positive. Each of the ten cluster-state comparisons has a coordinate difference excluding zero; uniqueness of the full row equations prevents later coalescence.

Equivariance supplies the remaining relative intervals exactly. Thus `T=6L=18 ell` is a positive return period. Only the three short segments are integrated for the proof; their matching conditions are themselves validated. No least-period assertion is needed for the stability conclusion.

## All twenty-three tangent directions

The exact invariant decomposition consists of the fourteen-dimensional quotient and one three-dimensional internal contrast for each of the three pairs. The singleton clusters contribute no internal contrasts. Hence `14+3*3=23` covers the whole relative-coordinate model.

Let H be the ordered product of the three quotient fundamental matrices. The relative return derivative is

```math
Q=P^{-1}H=P^5H,\qquad M_{\rm quotient}=Q^6.
```

The nonconstant autonomous orbit supplies an exact flow eigenvalue one for Q: equivariance gives `P^{-1}F(P y0)=F(y0)`. The characteristic polynomial is therefore divided by the exact factor `lambda-1`. An interval Schur test places every remaining root of Q inside radius 0.7. Spectral mapping puts the thirteen nonflow quotient multipliers for time T inside

```math
0.7^6=0.117649<0.13.
```

For pair j, let Bj be its internal three-dimensional contrast fundamental matrix over one relative interval L. Pair zero follows the cluster path `0,3,1,0,3,1`, so its full return is

```math
M_0=(B_1B_3B_0)^2.
```

The other two pairs use the corresponding cyclic products. All three full matrices are enclosed and tested separately, and their cubic characteristic polynomials pass the radius-0.13 Schur test. The fixed contrast decomposition then gives nine additional stable directions. Together with the quotient, this proves the stated one-neutral/twenty-two-stable count.

A floating-point real-Schur basis proposes a useful coordinate change for Q. Its entries are saved as exact dyadic rationals. Arb verifies invertibility and encloses the complete similarity transform, preserving every matrix radius. No approximate orthogonality assertion or numerical eigenvalue is used as proof. The Schur recursion and exact flow-factor argument are the same ones documented in the [four-pair certificate](motor-paired-certificate.md).

## Evidence, controls and limits

The [first successful certificate](../results/motor-five-cycle-proof/attempt-01.json) uses 448-bit Arb arithmetic, 144 subdivisions per shooting segment and Taylor order 48. The Krawczyk contraction is below `1.985e-19`; the largest inclusion ratio is below `1.078e-6`. Both the ordinary stability gates and the sharper bounds used above pass. The run took approximately 861 seconds, including three Newton center proposals. Newton supplies guesses only; existence and stability depend on the subsequent strict interval tests.

The [independent recomputation](../results/motor-five-cycle-proof/recheck-01.json) changes these settings to 512 bits, 160 subdivisions per segment and order 50. It performs no Newton steps: the earlier center is only a proposal, and the inclusion, complete spectrum and sharper bounds are recomputed. Every gate passes again, with contraction below `1.960e-19`; the run took approximately 511 seconds. The two shooting boxes have identical stored dyadic endpoints, and their shared executed sources and archived input hashes agree. The recheck uses the same interval method and arithmetic library, rather than a second proof implementation.

The [independent controls](../results/motor-five-cycle-proof/analytic-controls.json) check 6,240 weighted-field, Jacobian, full-space contrast and permutation components at rational test states; 2,844 optimized Taylor coefficients against the generic sparse recurrence; and 2,453 coefficients against an exact rational constant-matrix exponential oracle. They also verify the restriction's exact right inverse and the permutation's order six.

The new standalone [validator](../research/validate_motor_five_cycle.py) adapts the sparse field and multiple-shooting layout from the previously reviewed saddle implementation. No published proof module is modified or assigned different globals. The generic [similarity checker](../research/certify_motor_saddle_spectrum.py) is reused with the new fourteen-dimensional matrix. Executed source snapshots, exact input hashes and all integration diagnostics accompany each certificate.

Read-only review checked the weighted field, phase convention, multiple-shooting Jacobian, ordered products, relative-return identity, sharp spectral bounds and complete tangent-space count. The computational basis remains the custom Picard/Taylor integrator and Python-FLINT/Arb; this is not proof-assistant formalization or external peer review.

Choose unused output paths when reproducing:

```bash
python research/seed_motor_five_cycle.py --output outputs/five-seed
python research/validate_motor_five_cycle.py --seed outputs/five-seed/seed.json --output outputs/five-cycle.json
python research/check_motor_five_proof.py --output outputs/five-controls.json
python research/validate_motor_five_cycle.py --verify-input results/motor-five-cycle-proof/attempt-01.json --newton 0 --bits 512 --steps 160 --order 50 --output outputs/five-recheck.json
```

The seed generator uses the separate simulation dependencies, while the interval validator and algebra controls use `requirements-proof.txt`. The saved seed and dyadic basis can be used directly, so numerical seed generation need not be repeated to rerun the proof.

The [packet manifest](../results/motor-five-cycle-proof/manifest.json) lists the report, live sources, complete seed/control/proof evidence and their raw-byte SHA-256 hashes. It excludes itself from its own inventory. Both complete integrations retain all executed source snapshots.

## Symmetry copies and density lift

Distinct cluster states make the pointwise stabilizer exactly `H0=(C2)^3`, of order eight: only the two rows within each of the three pairs can be exchanged without moving the state. Any permutation preserving the whole orbit commutes with the flow and therefore acts by a phase translation. Its kernel is H0, so the finite orbital quotient is cyclic and lies in

```math
N_{S_8}(H_0)/H_0\simeq S_3\times S_2.
```

The largest cyclic subgroup of `S3 x S2` has order six, and the certified paired-cluster 3-cycle together with the singleton swap supplies an element of that order. The full orbital stabilizer therefore has order `8*6=48`, giving exactly

```math
\frac{8!}{48}=840
```

distinct labeled orbit copies. This argument uses the existence of a minimal period abstractly, not an assertion that the displayed T is minimal. These 840 cycles are symmetry copies of one attracting orbit class, not 840 new mechanisms or independent novelty claims.

An [independent integer enumeration](../results/motor-five-cycle-proof/isotropy-control-02.json), implemented in [check_motor_five_isotropy.py](../research/check_motor_five_isotropy.py), checks all 40,320 row permutations, the 96-element normalizer, quotient orders, the 48-element generated orbital group, and all 840 cosets. It references the successful interval certificate by hash. The earlier output is also retained; the second version corrects a left/right coset label without changing the groups, counts or argument.

For this fixed N=8, the [density-PDE lifting theorem](motor-pde-lifting.md) adds only the stable omitted-mode return factor `exp(-T)`, so the periodic density orbit inherits attraction under the stated sinusoidal-rate assumptions. Stability is not claimed after replication to more rows: duplicating the two singleton clusters would introduce additional internal contrast blocks that this eight-row certificate does not test.

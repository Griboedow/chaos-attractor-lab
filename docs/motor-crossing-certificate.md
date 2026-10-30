# A unique, nonzero-speed transverse crossing

The four-plus-four motor cycle has a **unique local transverse multiplier crossing with exactly six critical directions** in the narrow parameter interval certified below. Its crossing speed is strictly nonzero. This strengthens the earlier [periodic-orbit and intermediate-crossing certificate](motor-cycle-certificate.md), which established at least one crossing but did not certify uniqueness, speed or exact multiplicity there.

This remains a local theorem about a periodic orbit. It does not prove a nonlinear flip bifurcation, a branch of doubled cycles, a route to chaos or a chaotic attractor. The sixfold degeneracy excludes a direct application of the generic simple-multiplier flip theorem.

## Exact statement

The model, time normalization and exact algebraic gain are those in the [cycle certificate](motor-cycle-certificate.md):

```math
r=1+\frac{19/20}{\sqrt{32001/500000}},\qquad b=1.
```

The finite decimals defining the interval endpoints mean exact rational numbers:

```math
I_*=[3.370218399579184091349949673383546139330,\quad
     3.370218399579184091349969673383546139330].
```

For every rho in this interval, a locally unique half-period-symmetric periodic solution exists in the stored common shooting box. There is exactly one rho in its interior at which a transverse Floquet multiplier equals minus one. At that parameter:

- Its algebraic and geometric multiplicities in the full 23-dimensional system are exactly six.
- Sixteen multipliers lie strictly inside the unit circle; the remaining multiplier is the exact flow multiplier one.
- The derivative of the critical multiplier with respect to rho satisfies the outward bounds below.

```math
6.733318 < \frac{d\mu_*}{d\rho} < 6.733320.
```

The other two roots of each transverse block remain in the open intervals `(-0.045, -0.044)` and `(0.587, 0.590)` throughout the family. The critical root remains in `(-1.001, -0.999)` and is simple within each three-dimensional block. Three row contrasts in each of the two clusters supply six independent copies.

The lower endpoint has six unstable, one neutral and sixteen stable multipliers; the upper endpoint has zero unstable, one neutral and twenty-two stable multipliers. Increasing the source activity decreases rho, so activity crosses this boundary toward transverse instability.

The corresponding source activity and positive period have the enclosures

```math
\delta_* = 0.24147384540529444445082 \;\mathbin{\pm}\; 8.46\times10^{-24},
```

```math
T_* = 3.34200881521832174060 \;\mathbin{\pm}\; 5.31\times10^{-21}.
```

These are conservative readable enclosures; the stored dyadic midpoint/radius pairs are authoritative. Parameter digits describe the exact mathematical interpretation of the printed source constants, not experimental precision.

The new parameter and shooting boxes lie strictly inside the original validated family boxes. Consequently this is the same local periodic family. The independent [least-period certificate for that original family](../results/motor-period-proof/family.json) also applies here, making the displayed positive period its least period. Uniqueness of the crossing is established in `I_*`; this calculation does not exclude other crossings elsewhere in the wider original interval or on other periodic branches.

## How the derivative is enclosed

Let `q = (Y0, Z0, L)` denote the half-period shooting unknowns and let `F(q,rho)=0` be the symmetry closure equation. The implementation freshly validates a common shooting box with a uniform Krawczyk contraction. A parameter mean-value enclosure bounds the center residual using a validated partial rho sensitivity. It avoids repeatedly folding the same tiny parameter uncertainty into the state interval; no uncertainty is discarded.

Validated integration then encloses both the shooting Jacobian and its fixed-shooting partial parameter derivative. A rigorous ball-matrix solve gives

```math
q'(\rho)=-[D_qF]^{-1}\partial_\rho F.
```

For any augmented state or fundamental-matrix system written on normalized time as `x_u = L H(x,rho)`, total branch sensitivity satisfies

```math
z_u=L\left(D_xH\,z+\partial_\rho H+\frac{L'}{L}H\right).
```

The initial core sensitivity is `(0, Y0', Z0')`; initial fundamental-matrix sensitivities are zero. The code generates these degree-two equations by differentiating each sparse monomial. Coefficients are affine in rho, so their partial derivatives are exact coefficient differences, not differences of integrated trajectories. Every equation is integrated with the same Picard existence tubes and interval Taylor remainders used by the base certificate.

Half-period exchange gives the complete transverse monodromy `M_A = B_B B_A`; its derivative is formed by the product rule. The determinant function and its derivative are

```math
g(\rho)=\det(I+M_A(\rho)),\qquad
g'(\rho)=D\det(I+M_A)[M_A'].
```

The determinant derivative is evaluated by multilinearity, without inverting a matrix that becomes singular at the crossing. Uniform enclosures give `g' > 10.21745` for both transverse blocks. Newly integrated endpoint certificates give strictly opposite determinant signs. Continuity and strict monotonicity therefore establish exactly one zero in the open interval.

Opposite characteristic-polynomial signs at three disjoint interval pairs isolate all three real roots throughout the family. With `P(lambda,rho)=det(lambda I-M_A(rho))`, a separate enclosure of `-P_rho/P_lambda` is positive on the critical root interval. At the crossing it sharpens to the derivative bounds stated above. The two nonflow core roots satisfy strict Schur inequalities; the two passive cluster modes have modulus `exp(-T)`. This accounts for every one of the 23 multipliers.

## Evidence and reproduction

The successful [first certificate](../results/motor-crossing-proof/attempt-05.json) used 256-bit arithmetic, 160 subdivisions and Taylor order 26. The [recomputation](../results/motor-crossing-proof/recheck-01.json) used 320 bits, 192 subdivisions and order 28. The [final guarded recomputation](../results/motor-crossing-proof/recheck-zero-guard.json) repeats the latter settings after a coefficient-handling robustness correction. All pass their recorded gates. The uniform Krawczyk contraction is below `3.216e-17`. The first and second runs took approximately 77 and 103 seconds respectively on the local CPU.

The [analytic controls](../results/motor-crossing-proof/analytic-controls-zero-guard.json) check 408 field components against independently written core and transverse derivative formulas, a soluble sensitivity problem with a varying initial state and duration, a determinant derivative at a singular matrix, and retention of a nonzero-radius coefficient that contains zero. The maximum field-assembly discrepancy is below `1.7e-75` in ball arithmetic. These controls support software review; they do not replace the interval proof.

```bash
python -m pip install -r requirements-proof.txt
python research/validate_motor_crossing.py --bits 320 --steps 192 --order 28 --output outputs/crossing-recheck.json
python research/check_motor_crossing.py --output outputs/crossing-controls.json
```

Choose unused output paths. The scripts refuse to overwrite evidence. Each integration result records exact input balls, complete margins, endpoint records, source hashes and source snapshots. The approximate Newton refinements propose centers only; the final inclusion and derivative bounds are recomputed separately. The validator is [validate_motor_crossing.py](../research/validate_motor_crossing.py), with [independent controls](../research/check_motor_crossing.py).

All preliminary attempts remain in the [results directory](../results/motor-crossing-proof/). Attempts 01–03 failed because direct interval sensitivity propagation over the wider input box was too broad; attempt 02 additionally exposed nonfinite diagnostic serialization, which was corrected. Attempt 04 validated monotonicity on the narrow family but did not isolate endpoint signs in both blocks. Tightening only the endpoint shooting boxes produced attempt 05. None of these failures is evidence against existence of the cycle or crossing.

A read-only peer review checked the product-rule assembly, normalized-time duration term, implicit derivative solve, determinant derivative, common-box links and spectral isolation. It also found a general robustness issue in optional coefficient omission: an Arb inequality against zero must not replace an exact-zero predicate. The code now omits a coefficient only when `is_zero()` is true. In the earlier successful crossing runs, the duration coefficient was strictly negative and the coefficient derivatives were exact integers, so those particular calculations had not dropped an uncertain zero-containing coefficient. The final recomputation and added control verify the corrected implementation. The changed-discretization recomputation uses the same interval integrator; it is not an independent implementation, a proof-assistant formalization or an external peer review. The trusted computational base includes the custom Taylor/Picard code, Python-FLINT/Arb and their arithmetic/runtime stack.

## References

- [Original motor equations and exact normalization](motor-rows-source-audit.md); [base interval certificate and method references](motor-cycle-certificate.md).
- Python-FLINT developers, [real ball matrix operations](https://python-flint.readthedocs.io/en/latest/arb_mat.html). The default/`lu` solve encloses the result; the explicitly approximate solve is not used.
- I. Walawska and D. Wilczak, “An implicit algorithm for validated enclosures of the solutions to variational equations for ODEs,” *Applied Mathematics and Computation* **291**, 303–322 (2016). [Primary manuscript](https://arxiv.org/abs/1509.07388), [DOI](https://doi.org/10.1016/j.amc.2016.07.005). This is background on validated variational integration, not the implementation used here.

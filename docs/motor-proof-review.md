# Review of a proposed rigorous motor-cycle certificate

This is a proof design and independent review, not a completed interval certificate. The existing floating-point cycle and Floquet computations provide initial approximations. They do not by themselves establish the conclusions below with rigorous error bounds. The immediate target is the periodic four-plus-four cluster orbit at the source activity-path parameter, not a proof that the observed irregular set is a hyperchaotic attractor.

The proposed proof has three parts: validated integration, a phase-fixed shooting existence test, and certified multiplier counts. Each part must succeed before the corresponding theorem is claimed. The standard separation between a rough enclosure and a Taylor endpoint enclosure is described by Nedialkov, Jackson and Corliss [1]; rigorous variational integration is treated in [2]. The implementation-specific reductions and elementary spectral tests below are derived here for this model.

## Exact reductions before numerical proof

The three-state core is

```math
\begin{aligned}
h'&=[-\rho+r(1-Z)]h-rY,\\
Y'&=h(1-Z)-Y,\\
Z'&=hY-Z.
\end{aligned}
```

The five-state cluster restriction adds

```math
u'=-u+hv,\qquad v'=-v-hu+b.
```

For every fixed real $`b`$, the graph $`u=bY`$, $`v=b(1-Z)`$ is invariant. With $`a=u-bY`$ and $`d=v-b(1-Z)`$, direct subtraction gives

```math
a'=-a+hd,\qquad d'=-d-ha,\qquad
(a^2+d^2)'=-2(a^2+d^2).
```

Consequently every nonconstant periodic solution with positive period lies on this graph. Its two additional Floquet multipliers have modulus exactly $`e^{-T}`$. A rigorous three-state periodic-orbit certificate therefore lifts to an exact five-state cluster cycle, without a separate five-state shooting problem. This reduction concerns cluster cycles; it does not reduce the general 23-state dynamics.

The eight-row embedding and contrast blocks are those derived in [the transverse stability report](motor-cluster-transverse.md). They are algebraic identities. Their numerical spot checks are useful implementation checks, but the identities themselves justify the decomposition into a core block, two stable slave modes, three copies of each of two three-state transverse blocks.

## Declare which parameters the theorem concerns

An exact rational decimal such as $`r=4.755146047751339`$ and the value obtained from exact printed source constants are different mathematical parameters. The certificate must identify its choice. An interval enclosing the latter may be used for a uniform parameter-family certificate; it must not be silently replaced by its midpoint.

The source calibration of $`d`$ simplifies the activity-independent coefficient to

```math
r=1+\frac{k}{\eta c\sqrt{\Omega_0}}.
```

Thus, interpreting the printed constants as exact, one may construct

```math
\Omega_0=\frac{32001}{500000},\qquad
\frac{k}{\eta c}=\frac{19}{20},\qquad
\rho_{\rm target}=\frac{1291}{440},\qquad b=1.
```

This avoids cancellation through the original expressions containing $`\pi`$. It does not assign experimental certainty to the printed constants. A small enclosing parameter box is acceptable if every subsequent enclosure is uniform on that box. Uniform shooting inequalities then establish one nearby root for each fixed admissible parameter tuple, not a single root common to the whole box.

## A sufficient validated Taylor step

Let $`w'=F(w;p)`$ be the polynomial augmented system, with a fixed parameter tuple $`p`$ in a declared box $`P`$. The augmented state can contain the core and fundamental matrices. Let $`X`$ enclose all possible initial states for one step of exact positive length $`\Delta`$.

First find a finite box $`B`$ and certify

```math
X+[0,\Delta]F(B;P)\subset\mathrm{int}(B).
```

A conservative, explicit sufficient condition for the Picard contraction is also

```math
\Delta\,\sup_{w\in B,\,p\in P}\|D_wF(w;p)\|_W<1,
```

where a weighted maximum norm may reduce pessimism. These inequalities supply existence, uniqueness and a tube containing every solution throughout the step. Their failure is a failed validation attempt, not permission to use an unvalidated numerical trajectory as the tube. Parameters are fixed along each trajectory even when their values range over $`P`$.

Compute normalized solution Taylor coefficients by polynomial series arithmetic:

```math
a_0=w_0,\qquad
a_{j+1}=\frac{1}{j+1}[s^j]\,
F\left(\sum_{i=0}^{j}a_i s^i;p\right).
```

Interval evaluation at $`X`$ encloses the coefficients at all initial states. Evaluating the same coefficient recurrence through $`a_{n+1}`$ with $`a_0=B`$ bounds the normalized derivative appearing in Taylor's remainder at every possible intermediate point. Therefore

```math
\Phi_\Delta(X;P)\subset
\sum_{j=0}^{n} a_j(X;P)\Delta^j
+\Delta^{n+1}a_{n+1}(B;P).
```

All products, divisions, interval hulls and strict comparisons need directed enclosure arithmetic. Comparing order $`n`$ and order $`n+1`$ numerical answers is not a replacement for this remainder. Increasing arithmetic precision does not, by itself, bound truncation or wrapping error.

For a variable period, use the fixed interval $`s\in[0,1]`$ and $`dy/ds=T f(y)`$, with $`T`$ an uncertain positive parameter. A rational subdivision of this fixed interval avoids an uncertified floating-point terminal-time adjustment. The variational systems use the same scaled time and the same uncertain core trajectory.

## Wrapping and recentering

Direct interval Taylor propagation is valid if it completes, but can become too wide. A mean-value update must retain all initial uncertainty:

```math
\Phi_\Delta(X;P)\subset
\Phi_\Delta(c;P)+[D_y\Phi_\Delta(X;P)](X-c),\qquad c\in X.
```

Both the center image and derivative enclosure on the entire box are validated quantities. Resetting the state to a numerical center, resetting a fundamental matrix without accumulating its effect, or dropping a small previous radius destroys the proof.

A possible next implementation is a doubleton representation $`c+A\xi+R`$, where the original uncertainty $`\xi`$ is retained. Given a validated derivative enclosure $`[M]`$ on the hull of that set, exact selected centers $`c_+`$ and an exact selected matrix $`A_+`$, an enclosure for the new remainder is

```math
R_+=\bigl[\Phi_\Delta(c;P)-c_+\bigr]
+\bigl([M]A-A_+\bigr)\xi+[M]R.
```

Then $`c_++A_+\xi+R_+`$ encloses the next set. Any chosen coordinate transformation is admissible only if the transformed remainder and inverse are enclosed rigorously. Ordinary floating-point QR can propose coordinates; it cannot justify deleting their roundoff or uncertainty. Parameter dependence must also remain enclosed; augmenting constant parameters is one valid option.

If a single long shooting interval wraps excessively, multiple shooting is a clean alternative. With $`m`$ segments, use $`m`$ core states plus $`T`$, impose $`m`$ three-component matching equations and one phase equation, and enclose each shorter flow. This increases the algebraic system to $`3m+1`$ variables while shortening every validated IVP. It does not change the model.

## Phase-fixed shooting and Krawczyk inclusion

Use unknowns $`q=(Y_0,Z_0,T)`$ and initial state $`y_0(q)=(0,Y_0,Z_0)`$. Define

```math
G(q)=\Phi_T(y_0(q))-y_0(q).
```

Require $`T>0`$ and $`h'(0)=-rY_0>0`$ throughout the shooting box. The latter excludes an equilibrium and establishes a transverse phase section. It does not alone establish that $`T`$ is the least period; proving first return requires additional section-crossing exclusions. A positive period is sufficient for the immediate orbit and stability claims.

Validated integration of

```math
y_s=T f(y),\qquad M_s=T Df(y)M,\qquad M(0)=I
```

over all initial conditions and periods in the shooting box gives the needed derivative enclosure. Its columns are

```math
DG(q)=\left[
M e_Y-e_Y,\quad M e_Z-e_Z,\quad f(\Phi_T(y_0))
\right].
```

The last column is the derivative with respect to the original period, not the scaled-time field $`T f`$. A fundamental matrix computed only along the numerical center is insufficient to enclose $`DG`$ throughout the box.

Choose an exact center $`q_c`$, positive radii $`d_i`$, a box $`Q`$, and a nonsingular exact point preconditioner $`R`$. A floating-point inverse may suggest $`R`$; its stored dyadic entries can then be treated as exact. Form, with rigorous enclosures,

```math
K=q_c-R[G(q_c)]+\bigl(I-R[DG(Q)]\bigr)(Q-q_c).
```

An explicit sufficient certificate is strict containment $`K\subset\mathrm{int}(Q)`$ together with the weighted contraction bound

```math
\max_i\sum_j
\sup\left|\bigl(I-R[DG(Q)]\bigr)_{ij}\right|
\frac{d_j}{d_i}<1.
```

Banach's theorem then gives a unique root in this phase-fixed box. This is a conservative form of the standard verified nonlinear-system method [3]. Record strict inclusion and contraction margins, not just a success flag. All proof decisions must use enclosing bounds; converted Python floats are suitable for display only.

## Certified stability without a 23-state interval eigensolver

After the periodic solution is enclosed, enclose the monodromy matrices for every initial state and period in its certified box. Reusing a point-orbit monodromy would leave a gap.

For the core matrix $`M_c`$, exact periodic closure and a nonzero vector field imply an exact multiplier $`1`$. The other two roots satisfy

```math
z^2-sz+p=0,\qquad s=\mathrm{tr}(M_c)-1,\quad p=\det(M_c).
```

Strict interval inequalities

```math
1-p>0,\qquad 1-s+p>0,\qquad 1+s+p>0
```

certify that both roots lie inside the unit disk. For a conjugate pair this follows from $`p=|z|^2<1`$. For real roots, positivity of the polynomial at $`1`$ and $`-1`$, together with $`p<1`$, rules out every configuration with a root outside the unit interval. The exact neutral factor is used analytically; an interval eigenvalue merely overlapping $`1`$ is not declared equal to $`1`$.

For each real transverse matrix $`M_A,M_B`$, evaluate $`P(z)=\det(zI-M)`$ at rational endpoints of three disjoint real intervals suggested by the numerical data, for example

```math
[-3.71,-3.70],\qquad[-0.012,-0.010],\qquad[0.420,0.423].
```

If rigorous opposite endpoint signs hold for each interval, the intermediate value theorem gives at least one real root in each. A degree-three polynomial then has exactly these three roots, counted with multiplicity, and each is simple. One is outside the unit disk and two are inside. These intervals are proposed tests, not currently certified enclosures.

Perform this test independently for A and B. Exact zero-sum contrast identities repeat each block three times. The core contributes one neutral and two stable directions, and the slave modes contribute two stable directions. Successful tests therefore establish **six unstable, one neutral and sixteen stable directions of this periodic orbit**. Equal A/B spectra or an exact half-period exchange symmetry are not required for this count.

FLINT also supplies certified complex eigenvalue routines and multiplicity-aware variants [4]. Approximate eigenvalues, an arbitrary tolerance, or an enclosure covering the spectrum without assigning root counts are insufficient substitutes. Python bindings distinguish rigorous eigenvalue isolation from the explicitly approximate algorithm [5].

## What this certificate would not prove

It would not prove a strange attractor, asymptotic hyperchaos of the observed trajectory, a connecting orbit, a global route from cluster breakup, or scientific priority. An unstable periodic orbit with six expanding directions is a local dynamical result.

A rigorous transverse crossing is a separate task. One possible extended system combines core shooting with $`(M_A+I)v=0`$ and a normalization of $`v`$, adding the parameter $`\rho`$ as an unknown. Certifying its derivative requires validated parameter and second-order variational information or another rigorous derivative bound. Ordinary centered finite differences are not such a bound. Even a certified crossing with a nonzero slope would not settle the nonlinear, symmetry-degenerate flip bifurcation or establish a period-doubling cascade.

## Implementation and evidence requirements

Arb real balls include arithmetic rounding in their radii, but this does not supply an ODE truncation theorem automatically. Use integer/rational strings or explicitly declared dyadic values for inputs; decimal-to-binary conversions must preserve the intended enclosure. The documented upper/lower absolute bounds are appropriate for certified inequalities [6].

Retain exact parameter definitions, precision, package versions, source hashes, shooting box, time subdivision, Taylor order, Picard inclusion margins, remainder bounds, all failed validations, final root inclusion and spectral-sign margins. A separate verifier should recompute the inequalities. Test implementations against exactly soluble linear equations and the symbolic model identities, while clearly distinguishing such tests from a certificate for this nonlinear cycle.

No validated IVP, Krawczyk inclusion or spectral sign test was run as part of this review.

## References

1. N. S. Nedialkov, K. R. Jackson and G. F. Corliss, “Validated solutions of initial value problems for ordinary differential equations,” *Applied Mathematics and Computation* **105**, 21–68 (1999). [DOI and primary article](https://doi.org/10.1016/S0096-3003(98)10083-8).
2. I. Walawska and D. Wilczak, “An implicit algorithm for validated enclosures of the solutions to variational equations for ODEs,” *Applied Mathematics and Computation* **291**, 303–322 (2016). [Primary preprint](https://arxiv.org/abs/1509.07388); [DOI](https://doi.org/10.1016/j.amc.2016.07.005).
3. S. M. Rump, “Verification methods: rigorous results using floating-point arithmetic,” *Acta Numerica* **19**, 287–449 (2010). [Author-hosted paper](https://www.tuhh.de/ti3/rump/intlab/ActaNumerica2010.pdf); [author's explanation of the inclusion methods](https://www.tuhh.de/ti3/rump/Research_Rump/topics.shtml).
4. FLINT developers, [certified eigenvalue algorithms and their count guarantees](https://github.com/flintlib/flint/blob/main/doc/src/acb_mat/eigenvalues.md).
5. python-flint developers, [`acb_mat.eig` documentation](https://python-flint.readthedocs.io/en/latest/acb_mat.html).
6. python-flint developers, [real-ball construction and directed bounds](https://python-flint.readthedocs.io/en/latest/arb.html).

## Completed implementation review — 30 October 2026

The design above preceded the implementation. A subsequent read-only review of [the validator](../research/validate_motor_cycle.py) and its stored certificates found the final mathematical structure consistent with that design. The completed result and reproduction instructions belong in [the motor-cycle certificate report](motor-cycle-certificate.md). This review did not run another scientific integration.

The implementation uses half-period shooting instead of the full-period formulation above. With $`R=\mathrm{diag}(-1,-1,1)`$, the core is equivariant and the certified equation is $`\Phi_H(0,Y,Z)=R(0,Y,Z)`$. Thus the full period is $`2H`$. The shooting derivatives have the required $`+e_Y`$, $`-e_Z`$ and unscaled endpoint-field terms. Full monodromies are reconstructed as $`(R M_c)^2`$ for the core and $`M_B M_A`$, $`M_A M_B`$ for the transverse blocks. The jet recurrence, weighted Picard bound, parameter/duration enclosures and these compositions were inspected directly.

Two issues found during review were corrected: maxima of rounded Arb row sums now use directed upper bounds before comparison, and invalid negative Taylor orders and other invalid numerical arguments are rejected. Brackets touching $`-1`$ correctly classify their roots strictly inside their open interiors because the endpoint signs are strict. The fallback bracket straddling $`-1`$ does not certify a stable/unstable count.

The [target recomputation](../results/motor-proof/target-recheck.json) reports the complete 6 unstable / 1 neutral / 16 stable count. The [uniform crossing-family certificate](../results/motor-proof/crossing-family-attempt-01.json) has contraction bound below 0.032351. The [lower endpoint](../results/motor-proof/crossing-low.json) certifies 6 / 1 / 16, and the [upper endpoint](../results/motor-proof/crossing-high.json) certifies 0 / 1 / 22. Both transverse determinants $`\det(I+M)`$ are strictly negative at the lower endpoint and strictly positive at the upper endpoint. Exact rational comparisons of the stored dyadic endpoints independently confirmed that both endpoint shooting boxes lie strictly inside the original family box and both endpoint parameter enclosures lie inside its parameter enclosure. This is the branch-identification check needed for the continuity argument.

One serialization qualification matters: the JSON midpoint/radius pairs are lossless, but constructing an Arb ball from those pairs can round its stored radius upward. Such a reconstruction is enclosing, not necessarily bitwise identical. Recomputing on the enlarged box is conservative. Containment in a previously certified family box must instead use its exact stored dyadic endpoints, or another inward-safe comparison; testing only against a widened reconstruction of that family box would leave a gap. The independent branch-membership check above used exact rational endpoints.

The uniform family and opposite endpoint determinant signs establish existence of at least one parameter with a transverse multiplier $`-1`$. They do not establish a unique crossing, a nonzero crossing slope, a nonlinear flip bifurcation or a route to the observed hyperchaotic candidate. The local periodic-orbit certificates do not constitute a proof of a hyperchaotic attractor. An independent recomputation is a further reproducibility check; this code review is not a formal proof-assistant verification of the implementation or arithmetic library.

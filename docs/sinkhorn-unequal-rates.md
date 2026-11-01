# Unequal-rate matrix scaling: spatial numerical candidates

Two three-state matrix-scaling iterations produce separated curved bands.
They are numerical candidates of an explicitly modified Sinkhorn iteration:
the row and column updates use different fixed relaxation rates, with the
row rate above 2. A complete 4,864-input campaign retains the failures as well
as the selected cases. Chaotic attraction and historical uniqueness remain
unproved.

[Rotating and static portraits](../README.md#matrix-scaling-with-unequal-rates) ·
[Numerical packet and reproduction](../results/sinkhorn-unequal-rates/README.md) ·
[Earlier single-rate study](sinkhorn-screen.md)

## Map and input design

Let K be a strictly positive m × n matrix, L its entrywise logarithm, and
let a and b be its row and column sums. Thus u = v = 1 solves the scaling
problem exactly. We use the projective logarithms

```math
x_i=\log(u_i/u_0),\qquad y_j=\log(v_j/v_0),\qquad x_0=y_0=0.
```

There are m + n − 2 stored coordinates. Define

```math
\begin{aligned}
f_i(y)&=\log a_i-\log\sum_j e^{L_{ij}+y_j},&F_i(y)&=f_i(y)-f_0(y),\\
g_j(x)&=\log b_j-\log\sum_i e^{L_{ij}+x_i},&G_j(x)&=g_j(x)-g_0(x).
\end{aligned}
```

The update uses the new row state in its column step:

```math
x^{+}=r_xx+\omega_x F(y),\qquad
y^{+}=r_yy+\omega_y G(x^{+}),\qquad
r_x=1-\omega_x,\quad r_y=1-\omega_y.
```

Equal rates recover the published overrelaxed iteration of
[Lehmann et al.](https://arxiv.org/html/2012.12562v2). Unequal fixed rates
are an explicit extension here. We do not identify it with ordinary Sinkhorn
or with a safeguarded optimization method. The implementation uses stable
log-sum-exp and analytic tangent propagation, without clipping or noise.

The first unequal-rate screen uses independent normal log entries. Some
resulting row masses are exponentially small, and several spatial projections
are almost planar. The subsequent 3 × 2 screen keeps all row masses comparable:

```math
K_{i0}=c_i(1-s_i),\qquad K_{i1}=c_i s_i,\qquad
s_i=\frac{1}{1+e^{-R_i}},\qquad \sum_i c_i=1.
```

The three log odds R include a moderate value and values of opposite signs;
normalized lognormal row masses have log standard deviation 0.5. All exact
binary64 inputs, seeds and initial states are retained. Case numbers are
identifiers within a screen and matrix shape; Case 150 below means the
balanced-row 3 × 2 input, not an entry of another screen.

## Selected numerical candidates

Cases 150 and 987 have rates (2.1, 1.8) and (2.02, 1.9), respectively.
Their volume factors are −0.968 and −0.93636. Their fixed point at the origin
is unstable: the Jacobian has the unpaired eigenvalue 1 − row rate, whose
modulus exceeds one. Both have two positive finite-window QR rates.

Both trajectories complete step 106,040,000. Across the last twenty
five-million-step blocks, Case 150 has leading rates between 0.009391 and
0.009480, and Case 987 between 0.020770 and 0.020944. Maximum absolute
coordinates over that extension are 2.029 and 1.304, respectively.

| Case | Full QR rates after step 106,040,000 | Leading rate at 70 digits | Leading rate at 100 digits |
| --- | --- | ---: | ---: |
| 150 | 0.009444, 0.003727, −0.045694 | 0.009392 | 0.009579 |
| 987 | 0.020820, 0.010796, −0.097371 | 0.021210 | 0.020911 |

The QR spectra use 200,000 observations after a further 10,000-step tangent
warmup. Every 10,000-step spectrum block satisfies the exact logarithmic
volume identity within 8.33e-16. The independent probability-domain
precision runs restart from step 6,040,000, discard 2,000 steps and observe
200,000 steps. Binary64 kernels, rates and checkpoints are converted exactly
through integer ratios; marginals are recomputed at working precision.
These are finite-window estimates, not certified asymptotic exponents.

The portraits use all 250,000 consecutive states ending at step 106,040,000.
An orthogonal PCA rotation changes the viewing orientation only; no coordinate
is omitted or independently rescaled. The third principal direction accounts
for 9.12% and 6.51% of total variance. Discrete iterates are rendered as points,
without connecting segments. The animation rotates this fixed data through
600 frames in 24 seconds.

Each state control uses 20,000 warmup and 200,000 observed steps. The
32 nearby starts use independent offsets at scales 1e-6, 0.001, 0.03 and 0.3.
The 16 fresh starts use scales 1e-6, 0.1, 1 and 5. Eight parameter controls
change each relaxation rate by −0.001, 0 or +0.001, excluding the unchanged pair.

| Control | Case 150 | Case 987 |
| --- | --- | --- |
| Original start | Positive growth | Positive growth |
| 32 perturbed endpoints | All positive | All positive |
| 16 fresh starts | 5 positive; 11 cross the escape threshold | All positive |
| 8 nearby rate pairs | 7 positive; 1 negative | All positive |
| Ordinary Sinkhorn from the original start | Numerical zero after 41 steps | Numerical zero after 42 steps |

All four near-origin starts escape for Case 150. These finite tests do not
certify an open basin. For both cases, the above-2 row rate also prevents
reuse of the global absorbing-box argument from the single-rate study.
Long survival, positive growth and a reproducible picture are insufficient
to replace a trapping proof or exclude an exceptionally long transient.

## Complete campaign and rejected signals

| Screen | Inputs | Numerical captures | Completed pilot windows | Escape threshold | Positive flags |
| --- | ---: | ---: | ---: | ---: | ---: |
| Rectangular, common rate | 1,536 | 1,003 | 533 | 0 | 20 |
| Independent log entries, unequal rates | 2,304 | 1,314 | 361 | 629 | 12 |
| Balanced rows, unequal rates | 1,024 | 0 | 245 | 779 | 57 |

The rectangular screen covers 2 × 3, 3 × 2, 2 × 4 and 4 × 2 matrices,
384 inputs each. It observes 40,000 steps after 20,000 warmup steps and flags
growth above 1e-5. Every flag is continued to step 6,060,000. All seven
three-state flags lose positive growth. Across all 20 flags, 17 show
nontrivial numerical recurrences with lags 3, 5, 7, 9, 11 or 21; one is near
the origin; one remains slowly drifting. Only 2 × 4 Case 326 retains
substantial growth. Its four-state QR rates are approximately
0.30300, 0.01256, −0.02220 and −0.31341. Its 70- and 100-digit checks also
retain positive growth, but its diffuse projection is not selected as a portrait.

The initial unequal-rate screen observes 20,000 steps after 10,000 warmup
steps, across 2 × 2, 2 × 3 and 3 × 2 matrices. Every completed-window flag
above 0.005 is continued. Case 464 of the 3 × 2 screen crosses the escape
threshold during the first million-step extension. The other eleven retain
positive growth through step 6,030,000. Most sampled three-state shapes
are thin or nearly planar.

The balanced-row screen observes 30,000 steps after 10,000 warmup steps.
All 57 flags receive a million-step continuation. Four cross the escape
threshold, and all 53 survivors receive another five million steps. One
more escapes; Case 405 settles into a four-step numerical recurrence.
The remaining 51 retain positive growth at step 6,040,000. The two selected
portraits are chosen for their geometry, not for having the largest exponent.

The numerical capture threshold is a maximum absolute anchored coordinate
below 1e-11; the unequal-rate escape threshold is 100,000. Crossing that
threshold is a recorded finite event, not a proof of divergence to infinity.
Partial windows after capture or escape are excluded from rate comparisons.
The recurrence test compares 4,096 final states at all lags through 2,048,
with tolerance 1e-9 times one plus the window's largest absolute coordinate.
These are recurrence observations, not least-period certificates.

An input-generation defect in the initial unequal-rate driver is retained
and accounted for: its volume test preceded swapping the two rates on odd
cases. Consequently 76 of its 2,304 actual inputs have volume factor of
modulus at least one. All twelve flagged inputs contract volume. The
balanced-row driver tests the actual rates without swapping them; every
one of its 1,024 inputs satisfies a volume-factor modulus below 0.98.

## Jacobian structure and the meaning of local stability

The two triangular updates give the exact determinant

```math
\det DT(z)=r_x^{m-1}r_y^{n-1}.
```

If neither rate equals one, the map is a global diffeomorphism. Recover the
old y from the second update with the known new x, then recover the old x
from the first. Volume contraction alone does not imply boundedness.

A useful spectral identity holds at every state, not only at the solution.
Form the positive joint matrix using the old y and the newly updated x:

```math
H_{ij}=e^{x_i^{+}}K_{ij}e^{y_j},\qquad
h_i=\sum_j H_{ij},\quad k_j=\sum_i H_{ij},\qquad
S_{ij}=\frac{H_{ij}}{\sqrt{h_i k_j}}.
```

The derivatives of the row and column log sums are the two conditional
probabilities of this same H. Before removing constants, their products
are similar to $`SS^{\mathsf T}`$ or $`S^{\mathsf T}S`$.
The top singular value is one, simple because H is strictly positive;
all the others are strictly below one. Removing each constant direction
removes precisely this top pair. Write

```math
A=DF(y),\qquad B=DG(x^{+}),\qquad
DT=\begin{pmatrix}
r_x I&\omega_x A\\
\omega_y r_x B&r_y I+\omega_x\omega_y BA
\end{pmatrix}.
```

For each remaining squared singular value $`\mu=\sigma_j(S)^2`$,
the paired eigenvalues satisfy

```math
\lambda^2-t\lambda+r_xr_y=0,\qquad
t=r_x+r_y+\omega_x\omega_y\mu,\qquad 0\leq\mu\lt1.
```

There are also m − n unpaired eigenvalues r_x when m exceeds n, or
n − m eigenvalues r_y when n exceeds m. If both rates lie strictly between
zero and two, the real quadratic satisfies

```math
\begin{aligned}
1-t+r_xr_y&=\omega_x\omega_y(1-\mu)>0,\\
1+t+r_xr_y&=(1+r_x)(1+r_y)+\omega_x\omega_y\mu>0,\\
1-r_xr_y&>0.
\end{aligned}
```

The quadratic unit-disk criterion therefore puts every instantaneous
Jacobian eigenvalue strictly inside the unit disk. This is a spectral
condition, not contraction in a common norm. Products of different
Jacobians can still amplify perturbations. The statement applies to
the common-rate rectangular screen; it does **not** apply to the selected
above-2 unequal-rate portraits.

The factorization was checked independently at 240 moderate random states
of eight matrix shapes. Maximum eigenvalue matching error was 9.59e-13;
maximum determinant error was 8.89e-16. Separate checks at 18 actual retained
states compare the implementation with 100-digit probability-domain
arithmetic, centered finite differences and inverse reconstruction.
Maximum state and tangent errors are 2.77e-14 and 5.84e-15, respectively;
the finite-difference error is 1.17e-9. These checks support the
implementation; the universal spectral statement follows from the algebra.

## Prior art and remaining equivalence questions

Overrelaxation and its possible nonconvergence are established prior work.
[Thibault et al., §4](https://arxiv.org/html/1711.01851v2) also allow distinct
adaptive functions for the two projection steps, selected to decrease a
Lyapunov function. Our fixed above-2 rate does not satisfy their safeguard.
Two fixed block rates are already standard in
[modified SOR](https://www.sciencedirect.com/science/article/abs/pii/S0096300312010879).
Assigning different rates to two blocks is not itself a new mechanism.

The common-rate Jacobian at the solution has the block-SOR interpretation
in [Lehmann et al., §2](https://arxiv.org/html/2012.12562v2). The calculation
above uses the compatible joint matrix to extend the spectral identity to
every state. No priority claim is made for this elementary extension.
Nor is failure of global convergence despite stable instantaneous
eigenvalues new in general: explicit discrete Markus–Yamabe counterexamples
are given by [van den Essen and Hubbers, 1997](https://www.cs.ru.nl/E.Hubbers/pubs/Chaotic-Polynomial-Automorphisms-Counterexamples_1997_Advances-in-Applied-M.pdf).
Their displayed counterexamples involve escaping orbits, not these portraits.

Searches on 2026-11-01 combined Sinkhorn, matrix scaling, unequal relaxation,
modified SOR, chaos, hyperchaos and chaotic attractor. No matching
classification of the two selected matrix-scaling configurations was
identified in the inspected primary material. This bounded search does not
establish historical uniqueness. The smooth invertible map is not excluded
from equivalence to a known invertible chaotic map; curved bands and two
positive finite rates cannot decide that question. A trapping certificate,
an invariant-set analysis and a stronger equivalence audit remain necessary.

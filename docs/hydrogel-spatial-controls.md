# Hydrogel filament: spatial dynamics and mesh artifacts

## Scope

The planar implementation has been extended to a spatial, inextensible rod
with two evolving spontaneous curvatures and quasistatic twist. The solver
recovers planar forces and second-order temporal convergence, and produces
spatial beating and hook-shaped rotation at the example parameters of
[Boiardi et al., arXiv:2604.09486v1](https://arxiv.org/html/2604.09486v1).
The computations are independent segment approximations, not a reproduction
of the authors' COMSOL discretization.

The controls expose a substantial mesh artifact: a 16-segment approximation
admits a bent equilibrium absent from the continuum model. Coarse irregular
trajectories also encounter nonlinear-solver failures. These outcomes prevent
promotion of the coarse pictures as new attractors. A stronger-field spatial
scout completes but develops distances between separated rod segments below
the nominal filament diameter. Neither chaos nor historical priority is
established.

[Numerical packet and reproduction](../results/hydrogel-spatial-controls/README.md) ·
[Earlier planar controls](hydrogel-planar-controls.md).

The [longer refinement and variational study](hydrogel-variational-checks.md)
tests the thinner circular regime, checks perturbation growth and shows why
positive finite growth during slow drift does not establish a chaotic attractor.

## Spatial discretization

The source retains translational resistive-force drag, neglects inertia and
does not add a rotational drag law. Its moment is proportional to the
difference between visible and spontaneous strain. The two spontaneous
curvatures obey source equation (25); spontaneous twist is zero. The
[planar report](hydrogel-planar-controls.md#model-and-source-correspondence)
defines the common dimensionless parameters and drag coefficients. The
remaining torsional coefficient is

```math
\beta_3=\frac{\pi\eta^3}{2\lambda^4(1+\eta^2)(1+\nu)}.
```

Segment $`i`$ has orientation $`R_i\in\mathrm{SO}(3)`$ and length $`h=1/N`$.
Its tangent is the third column of $`R_i`$. Positions are accumulated from
these unit tangents, enforcing the discrete lengths exactly. Unit quaternions
represent the rotations. The clamp orientation is the identity; the first
bending interval has length $`h/2`$ and subsequent intervals length $`h`$.
As in the planar approximation, the first straight segment only approximates
the pointwise clamped tangent through its half-cell bending energy.

For each relative rotation, let $`v_i`$ be its principal rotation vector and
$`\ell_i`$ its bending-interval length. The midpoint director frame is the
geodesic midpoint of the neighboring frames. With $`\kappa_{i,3}=0`$,

```math
u_i=\frac{v_i}{\ell_i},\qquad
E_N=\frac12\sum_i\ell_i
(u_i-\kappa_i)^T
\mathrm{diag}(\beta_1,\beta_2,\beta_3)(u_i-\kappa_i).
```

The mechanical variation holds the internal curvatures fixed. Their evolution
uses the midpoint directors:

```math
\dot\kappa_1=-\kappa_1-\chi_1(d_2)_z,\qquad
\dot\kappa_2=-\eta^{-2}\kappa_2+
\eta^{-3}\chi_1(d_1)_z.
```

The local translational resistance is

```math
D_i=\xi R_i\,\mathrm{diag}(\zeta_1,\zeta_2,\zeta_3)R_i^T.
```

The discrete mechanical residual is the elastic orientation gradient plus
the drag virtual work. Drag is integrated exactly along each straight segment.
It has no component for rotation about that segment's tangent. Thus the axial
moment equation remains algebraic; no artificial rotational damping is used
to make an ordinary differential equation.

Constant-step BDF2 is applied to positions and spontaneous curvatures, with
one backward-Euler startup step. The mechanical equations are solved at the
new configuration. The linear curvature update is eliminated analytically,
but its dependence on the new directors is included in Newton's derivative.
Initial twist is adjusted to satisfy axial moment equilibrium while keeping
the initial centerline fixed.

The initial rotation vector at each segment center is

```math
\bigl(-0.01s(2-s),\;0.007s^2(3-2s),\;0\bigr),
```

with zero spontaneous curvature; planar controls omit its second component.
This is our specified perturbation, not a reconstruction of the paper's
initial velocity field. The principal logarithm requires resolved relative
rotations. Approaching its angle-$`\pi`$ branch or losing Newton convergence
invalidates a completed-trajectory claim.

## Algebra and time integration checks

The original implementation forms a complete central-difference Newton
matrix. The faster implementation uses the analytic drag derivative and
local differences of the two elastic hinges affected by each rotation.
Both sources are retained. A direct comparison gives:

| Check | Recorded discrepancy |
| --- | ---: |
| Elastic gradient versus central differences of the energy | `5.82e-11` relative |
| Drag virtual work versus independent three-point quadrature | `4.42e-16` relative |
| Spatial versus earlier planar elastic force | `5.37e-15` relative |
| Spatial versus earlier planar curvature relaxation | `1.78e-15` absolute |
| Analytic/local Newton matrix versus full differences, BDF1 | `4.38e-11` relative |
| Analytic/local Newton matrix versus full differences, BDF2 | `5.58e-11` relative |
| One converged step, original versus faster implementation | `1.12e-16` maximum state difference |

These are numerical implementation checks. They do not constitute an
independent continuum solver or an interval proof.

At 16 segments, restricting the spatial solver to the $`y,z`$ plane reproduces
the earlier periodic control at $`\chi_1=100`$, $`\lambda=50`$,
$`\eta=0.9`$, $`\xi=5\times10^{-6}`$ and $`\nu=0.5`$.

| Time step | Maximum angle difference from adaptive planar reference over times 0–20 | Maximum curvature difference | Mean period after time 15 |
| --- | ---: | ---: | ---: |
| 0.002 | 0.025893 | 0.196432 | 0.735924344 |
| 0.001 | 0.006500 | 0.049303 | 0.735866301 |

Halving the step reduces both history errors by a factor of about 3.98,
consistent with second-order convergence. The adaptive reference's asymptotic
period is 0.735846775. Its saved states are interpolated cubically; this
comparison is not a comparison of exact solutions. Phase drift makes the
maximum history discrepancy more sensitive than the final orbit's amplitude.

## A bent equilibrium created by the mesh

For the continuum source model under a constant axial field, every equilibrium
is straight. Indeed, zero velocity and a free end imply zero contact force
everywhere. Moment balance then gives zero moment everywhere, hence
$`u=\kappa`$. At equilibrium,

```math
\kappa_1=-\chi_1(d_2)_z,\qquad
\kappa_2=\eta^{-1}\chi_1(d_1)_z,\qquad
\kappa_3=0.
```

Substitution into the director kinematics gives a smooth first-order initial
value problem in arc length. Its clamped identity frame has the constant
solution, and uniqueness excludes any other equilibrium. This argument
assumes the source's straight reference rod, fixed axial field and stated
free-end conditions; it does not apply to naturally curved or externally
loaded rods.

The midpoint segment scheme can violate that property. Consider a circular
cross-section, a planar configuration and $`\chi_1=-\alpha<0`$. Vanishing
discrete moment requires

```math
\theta_0=\frac{\alpha}{2N}\sin\frac{\theta_0}{2},\qquad
\theta_i-\theta_{i-1}
=\frac{\alpha}{N}\sin\frac{\theta_i+\theta_{i-1}}2.
```

If $`2N<\alpha<\pi N`$, the first equation has only the zero solution, while the
first interior equation admits nonzero roots. Its derivative at zero is
$`1-\alpha/(2N)<0`$ and its value at $`\pi`$ is positive for the recorded
case. Subsequent angles can be obtained by continuous scalar root solves.
Each increment can be chosen in $`[-\alpha/N,\alpha/N]`$, inside the
principal-rotation interval under the stated inequality.
For $`N=16`$ and $`\alpha=33`$, this constructs an explicit bent discrete
equilibrium with tip approximately

```math
(0,\;0.06553205,\;-0.76875036).
```

The recorded scaled mechanical residual is below `5.51e-15` and the curvature
relaxation residual below `7.80e-14`. The 16-segment spatial run approaches a
rotated form of this stationary shape. In contrast, the 32- and 64-segment
runs retain hook-shaped rotation. The construction proves an unwanted
discrete equilibrium; it is not a new equilibrium of the physical model.

## Spatial controls

For the positive-field example, all three finer runs complete time 60 and
approach a spatial cycle. The table uses the last five upward tip-$`x`$
crossing intervals and the largest orientation return over the last ten
crossings. Relative rotations identify quaternion signs.

| Segments | Time step | Late mean period | Maximum orientation return angle |
| --- | ---: | ---: | ---: |
| 32 | 0.002 | 0.633728080 | `2.09e-5` |
| 32 | 0.001 | 0.633665837 | `2.17e-5` |
| 64 | 0.002 | 0.633637990 | `8.46e-6` |

The period changes by about 0.010% on halving the step at 32 segments and by
0.014% on doubling the mesh at step 0.002. These are convergence observations,
not rigorous error bounds. The corresponding 16-segment run fails Newton
convergence at time 41.742; its preceding irregular trace is not retained as
a resolved trajectory of the continuum system.

For the circular-section, reverse-field example, 32 and 64 segments give
late mean tip radii 0.08230 and 0.08441. Linear fits of the unwrapped azimuth
over times 80–100 give periods approximately 17.20 and 14.64, respectively.
The roughly 15% difference shows that the rotation's period is still poorly
resolved in space. These runs recover the qualitative rotating mode, not an
accurate continuum rotation rate. The coarse near-stationary run is separately
explained by the explicit discrete equilibrium above.

## Interpretation

At field $`\chi_1=150`$ with $`\lambda=50`$, $`\eta=0.9`$ and 64 segments,
the BDF2 scout completes times 0–60 with step 0.004. Its relative hinge angles
remain below 0.663 radians. Nonetheless, 67 of 417 sampled configurations over
times 40–60 have tested centerline distances below the nominal minor diameter
$`2\eta/\lambda=0.036`$. The minimum is approximately 0.000740 at time 41.968,
between segments 29 and 45 (zero-based indices). No contact force or nonlocal
hydrodynamic interaction is included in the source model. This trace is not
retained as a physically validated chaotic-attractor candidate.

The distance screen minimizes Euclidean distance between straight segments
whose entire material intervals are separated by at least $`4/\lambda`$.
It inspects configurations about 0.05 time units apart and omits closer
material neighbors. It is a finite geometry check, not a continuous-time
contact certificate or a finite-cross-section simulation. The corresponding
32-segment field-100 control has minimum tested clearance 0.08715 over its
late sampled window, above the same diameter threshold.

The positive-field periodic control and negative-field rotating control
recover the qualitative motions in the source's Figures 5a and 6c. Mesh
refinement is essential, particularly where curvature concentrates near the
clamp. Reported rotation periods and geometric statistics are numerical
observations; the authors' original trajectories were not available for
pointwise comparison.

## Circular-section symmetry and thinner filaments

A separate [three-run screen](../results/hydrogel-slender-screen/README.md)
uses a thinner circular rod: $`\lambda=127`$, $`\eta=1`$,
$`\xi=5\times10^{-6}`$ and $`\nu=0.5`$. Field 117 is the parameter set of
the source's Figure 6a and Supplementary Video 2; fields 130 and 150 increase
the drive. All three runs use 64 segments, step 0.006 and the same initial
perturbation as above. Their requested duration is 96. This is a small
declared screen, not a parameter-space exclusion.

### Removing a common rotation from state comparisons

For a circular section, both transverse bending stiffnesses and both
transverse drag coefficients coincide. With $`Z_\psi`$ denoting rotation
about the field axis, the equations and clamped identity frame are invariant
under

```math
R_i\mapsto Z_\psi R_i Z_\psi^T,\qquad
\kappa_i\mapsto Z_\psi\kappa_i,\qquad
r_i\mapsto Z_\psi r_i.
```

Here the two-component curvature is embedded as
$`\kappa=(\kappa_1,\kappa_2,0)`$. A changing laboratory-frame trace can
therefore reflect a common rotation. Define frame and curvature differences

```math
\Delta R_i(\psi)=R_i^A-Z_\psi R_i^B Z_\psi^T,\qquad
\Delta\kappa_i(\psi)=\kappa_i^A-Z_\psi\kappa_i^B.
```

We compare all director frames and curvatures by minimizing

```math
D^2(A,B)=\min_\psi\frac1N\sum_i
\left[\frac12\mathrm{tr}(\Delta R_i^T\Delta R_i)
+\frac{\Delta\kappa_i^T\Delta\kappa_i}{\chi_1^2}\right].
```

All differences inside the minimum are evaluated at the same angle $`\psi`$.
The trace term is half the squared Frobenius norm of the frame difference.

The angular part contains only first and second harmonics. Maximizing
$`a\cos\psi+b\sin\psi+c\cos2\psi+d\sin2\psi`$ reduces its stationary
angles to the real roots of

```math
(-b+2d)t^4+(-2a+8c)t^3-12dt^2+(-2a-8c)t+b+2d=0,
\qquad t=\tan(\psi/2).
```

The implementation evaluates every computed real root and the angles zero
and $`\pi`$, then chooses the largest score. This is an algebraic global
minimization reduced to a numerical polynomial solve, not a certified root
enclosure. Rotation matrices remove quaternion-sign ambiguity. The reduction
applies only to the circular section; it is not used for $`\eta=0.9`$.

A random-state residual check has relative equivariance error `2.55e-14`.
An imposed rotation of 0.731 radians is recovered with state distance below
`7.61e-17`; quaternion sign changes leave the comparison unchanged. A second
random pair agrees with a 4,097-angle grid to its sampling error. These are
implementation controls, not dynamical evidence.

Applied to the earlier 64-segment reverse-field rotation over times 80–100,
21 equally spaced states differ from the final state by at most 0.00800 after
alignment. For its one observed tip-height return, the distance falls from
1.89975 without alignment to 0.005344 with alignment. The small residual
shape variation and earlier mesh sensitivity remain; the calculation does
not certify an exact relative equilibrium or its rotation period.

### Three-run screen

The clearance screen examines 417 configurations in each run's final
20 saved time units, 0.048 units apart. It tests segment-index separations
of at least four, hence material gaps of at least $`3/64`$. The nominal
diameter is $`2/127\approx0.0157480`$.

| Field | Last saved time | Integration | Minimum tested clearance | Samples below diameter |
| ---: | ---: | --- | ---: | ---: |
| 117 | 96 | Completed | 0.0440251 | 0 / 417 |
| 130 | 96 | Completed | 0.000196665 | 137 / 417 |
| 150 | 52.608 | Newton failure at 52.620 | 0.000005531 | 315 / 417 |

For field 130 the minimum occurs at time 76.968 between segments 35 and 54;
for field 150 it occurs at time 37.176 between segments 4 and 49. Indices are
zero-based. Increasing slenderness has not removed close approaches from
these two stronger-field traces. Their irregularity is not promoted as a
physically validated chaotic-attractor candidate. The failed integration
also supplies no evidence about the model's later trajectory.

Full-state returns use upward crossings of the late mean tip height, an
observable unchanged by axial rotation. Quaternions are made continuous in
sign before cubic interpolation and renormalized at each event. Curvatures
are interpolated cubically. The record retains ordinary and aligned
distances at every tested lag from one to eight, where enough events exist.

The field-117 run has 21 such crossings in times 76–96. Its smallest
maximum aligned return distance among the eight lags is 0.8173, at lag three;
the final pair at that lag has distance 0.1024. Thus this finite window has
not established the published periodic regime. Slow drift of the mean height,
remaining transients and discretization
error can all affect these returns. They do not establish chaos or refute
the source's periodic example. This thinner-rod control requires longer
integration and separate temporal and spatial refinement before nearby
irregular traces can support a discovery claim.

## Prior art and remaining validation

The source already describes spatial periodic motion and its secondary
instability. Those features are not discoveries here. Searches for the exact
preprint identifier, DOI and hydrogel/morphoelastic chaos terms did not locate
an exact-model chaos analysis. Some matches concern chaotic printing of
hydrogels, which is a different mechanism. This bounded query set does not
establish absence of prior work, and broader active-filament chaos is known.

The remaining discovery task requires a persistent irregular regime with
positive tangent growth, temporal and spatial refinement, physically
interpretable geometry, and a model-specific priority/equivalence audit.
Neither a coarse irregular trace nor successful Newton solves satisfy those
requirements.

[Research index](../README.md#research-index) ·
[Candidate queue](candidate-queue.md)

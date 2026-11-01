# Hydrogel filament: planar solver controls

## Result and scope

An independent solver reproduces the two linear flutter thresholds of the
three-dimensional morphoelastic filament model of
[Boiardi et al. (2026)](https://arxiv.org/html/2604.09486v1).
A nonlinear planar control approaches a numerical periodic orbit under spatial
refinement and two implicit time integrators. An eight-case planar screen has
not established a new chaotic attractor: six sampled trajectories cross
themselves, and the other two approach numerical periodic orbits in longer,
refined checks.

This is a validation step toward the full spatial model. The computations
restrict motion to a principal plane and do not test stability to out-of-plane
perturbations. They introduce neither inertia nor rotational drag. No new
portrait or discovery-priority claim is attached to this packet.

[Code and complete numerical records](../results/hydrogel-planar-controls/README.md).

## Model and source correspondence

Arc length and time use the source's dimensionless variables. Write the tangent
and normal in either principal plane as

```math
\mathbf t=(\sin\theta,\cos\theta),\qquad
\mathbf n=(\cos\theta,-\sin\theta),\qquad
\mathbf r(s)=\int_0^s\mathbf t(u)\,du.
```

The clamp fixes position and angle at zero; the other end is force- and
moment-free. With effective spontaneous curvature $`k`$, bending energy and
translational resistance are

```math
E=\frac{\beta}{2}\int_0^1(\theta_s-k)^2\,ds,
\qquad
\mathbf D=\xi\bigl(\zeta_n\mathbf n\mathbf n^T+
\zeta_3\mathbf t\mathbf t^T\bigr).
```

Force balance is the weak overdamped equation

```math
\int_0^1\delta\mathbf r^T\mathbf D\,\dot{\mathbf r}\,ds
=-\beta\int_0^1(\theta_s-k)\,\delta\theta_s\,ds,
\qquad
\dot k=-a k-c\sin\theta.
```

The mechanical variation holds the internal variable $`k`$ fixed. The two
source planes correspond to the following coefficients; the sign change in
the first row follows from $`u_1=-\theta_s`$.

| Plane | Effective curvature | Bending coefficient | Normal resistance | Relaxation $`a`$ | Drive $`c`$ |
| --- | --- | --- | --- | --- | --- |
| $`y,z`$ | $`k=-\kappa_1`$ | $`\beta_1`$ | $`\zeta_2`$ | 1 | $`\chi_1`$ |
| $`x,z`$ | $`k=\kappa_2`$ | $`\beta_2`$ | $`\zeta_1`$ | $`\eta^{-2}`$ | $`\eta^{-3}\chi_1`$ |

The source definitions, implemented without fitted coefficients, are

```math
L=\log\frac{2\lambda}{1+\eta},\qquad
\zeta_1=\frac{4\pi}{L+\eta/(1+\eta)},\quad
\zeta_2=\frac{4\pi}{L+1/(1+\eta)},\quad
\zeta_3=\frac{2\pi}{L-1/2},
```

```math
\beta_1=\frac{\pi\eta}{4\lambda^4},\qquad
\beta_2=\frac{\pi\eta^3}{4\lambda^4}.
```

These relations follow from source equations (22), (23), (25), (36) and (39).
The full three-dimensional theory additionally has twist and director
constraints. Its planar invariant subspaces have the older morphoelastic
mechanism studied by
[Boiardi and Marchello (2024)](https://link.springer.com/article/10.1007/s11012-024-01765-7).
That paper assumes normal-to-tangential drag ratio two; the finite-slenderness
coefficients above retain the 2026 formulation. A planar simulation therefore
does not by itself demonstrate a new physical mechanism.

## Independent linear checks

Linearizing transverse displacement $`y`$ gives

```math
\xi\zeta_n\dot y=-\beta(y_{ssss}-k_{ss}),\qquad
\dot k=-a k-cy_s.
```

Clamped displacement and slope are essential boundary conditions; zero moment
and shear are natural conditions in the weak formulation. Cubic Hermite beam
elements give mass, bending and mixed matrices

```math
M_{ij}=\int\phi_i\phi_j\,ds,\quad
K_{ij}=\frac{\beta}{\xi\zeta_n}\int\phi_i''\phi_j''\,ds,\quad
H_{ij}=\frac{\beta}{\xi\zeta_n}\int\phi_i''\phi_j'\,ds.
```

Eliminating curvature gives the quadratic eigenvalue problem

```math
\bigl[\omega^2M+\omega(aM+K)+aK+cH\bigr]v=0.
```

This spatial discretization is independent of the angle-segment method used
for nonlinear evolution. At $`\lambda=50`$, $`\eta=0.9`$ and
$`\xi=5\times10^{-6}`$, the positive-field flutter thresholds are:

| Elements / segments | Hermite, $`x,z`$ | Angle, $`x,z`$ | Hermite, $`y,z`$ | Angle, $`y,z`$ |
| --- | ---: | ---: | ---: | ---: |
| 8 | 38.52849 | 42.20154 | 50.70481 | 56.08058 |
| 16 | 38.43688 | 39.35836 | 50.55510 | 51.89808 |
| 32 | 38.43100 | 38.66640 | 50.54536 | 50.88722 |
| 64 | 38.43062 | 38.49038 | 50.54474 | 50.63136 |

The article gives approximately 38 and 51 in Section 4.1. This comparison
recovers the stated thresholds, not the authors' unpublished mesh data.
For zero drive, the first passive cantilever eigenvalue is
$`-0.09990443390`$ analytically and $`-0.09990443560`$ with 32 Hermite elements.
An independent quadrature assembly agrees with the nonlinear drag matrix to
$`6.36\times10^{-16}`$ relative error. The nonlinear right-hand side agrees
with its separately assembled linearization to $`1.93\times10^{-12}`$ in the
recorded directional check. These are floating-point controls, not interval
certificates.

## Nonlinear discretization and periodic control

The nonlinear approximation has $`N`$ straight segments of length $`h=1/N`$,
with tangent angles at segment centers. The clamp-to-first-center bending
interval has length $`h/2`$; subsequent intervals have length $`h`$. If $`B`$
is the angle-difference matrix and $`W`$ contains those interval lengths,

```math
E_N=\frac{\beta}{2}(B\theta-k)^TW(B\theta-k),\qquad
M(\theta)\dot\theta=-\beta B^TW(B\theta-k).
```

Curvature relaxation samples the midpoint angle of each bending interval.
The first straight segment approximates the clamped tangent through its
half-cell bending energy; its pointwise boundary angle is not exactly zero at
finite mesh size. The prescribed continuum angle is recovered under refinement.
The drag mass integrates the position Jacobian exactly on each straight
segment. In particular, the current-segment factors integrate to $`1/2`$ and
$`1/3`$; replacing them by a midpoint lump would change the method. Refinement
is necessary: a finite segment chain is an approximation to the continuum rod.

For the source Figure 4a parameters $`\chi_1=100`$, $`\lambda=50`$,
$`\eta=0.9`$, $`\xi=5\times10^{-6}`$, motion is restricted to the $`y,z`$
plane. All runs start with $`\theta(s)=0.01s(2-s)`$ and $`k=0`$, which is our
declared perturbation rather than a reconstruction of the source's initial
velocity field. Statistics use times 30–60. Upward zero crossings of the tip's
transverse coordinate are located by each solver's event interpolant.

| Segments | Integrator, relative tolerance | Mean period | Largest full-state return error in last ten returns |
| --- | --- | ---: | ---: |
| 16 | BDF, `1e-7` | 0.735846775 | `2.53e-6` |
| 32 | BDF, `1e-8` | 0.734242007 | `2.32e-7` |
| 32 | Radau, `1e-8` | 0.734242009 | `5.72e-9` |
| 64 | BDF, `1e-8` | 0.733944369 | `2.18e-7` |

At a fixed mesh, the periods from BDF and Radau differ by $`1.59\times10^{-9}`$.
The 32-to-64 mesh period change is approximately 0.041%; the maximum sampled
transverse tip displacement changes from 0.11342 to 0.11327. No proper
centerline crossing was found in the 1,501 sampled late configurations of each
control. This establishes a numerical periodic control within the planar
restriction; it is not a proof of a three-dimensional attracting cycle.

## Finite exploratory slice

The eight declared cases use 32 segments, BDF with relative tolerance `1e-6`,
duration 80 and the same small initial bend. The table reports proper
intersections between nonadjacent straight segments in 2,001 configurations
sampled at interval 0.02 over times 40–80. Tangencies, finite thickness and
contact between samples are not covered by this test.

| Case | Slenderness $`\lambda`$ | Aspect ratio $`\eta`$ | Field $`\chi_1`$ | Configurations with crossings |
| --- | ---: | ---: | ---: | ---: |
| 0 | 50 | 0.9 | −33 | 0 / 2,001 |
| 1 | 50 | 0.9 | 150 | 747 / 2,001 |
| 2 | 50 | 0.9 | 300 | 1,481 / 2,001 |
| 3 | 127 | 1 | 117 | 0 / 2,001 |
| 4 | 127 | 1 | 200 | 1,937 / 2,001 |
| 5 | 127 | 1 | 400 | 1,998 / 2,001 |
| 6 | 200 | 1 | 160 | 377 / 2,001 |
| 7 | 200 | 1 | 300 | 1,994 / 2,001 |

Case numbers identify rows of this finite survey. All eight integrations
completed. These are planar trajectories, including where the source presents
three-dimensional motion at the same parameters. The six crossing trajectories
are not retained as physically admissible attractor candidates; they do not
exclude a nonintersecting state elsewhere in parameter space or in the full
spatial model. No Lyapunov exponent was inferred from return irregularity.

The short Case 0 run has a period standard deviation of 0.0659 and a large
return discrepancy. Extending the 32-segment run to time 600 with relative
tolerance `1e-8` gives period 2.414326460 and a last-ten return discrepancy below
`5.35e-7`. Its initial irregularity is therefore rejected as evidence of a
persistent chaotic attractor. A 64-segment Radau run to time 600 also closes,
with period 2.396433827 and return discrepancy below `2.00e-8`. The period still
changes by about 0.75% between these meshes; periodic closure is much better
resolved than the continuum value of the period.

For Case 3, a 64-segment BDF run to time 400 gives period 1.410320380 and a
last-ten return discrepancy below `1.62e-6`. No proper crossing was detected in
the late 15,001 samples of either Case 0 extension or the 10,001 late Case 3
samples. These finite checks reject the two short irregular histories as
retained chaos candidates, without proving global exclusion of chaos.

## Priority and remaining work

The [spatial follow-up](hydrogel-spatial-controls.md) now implements twist
equilibrium, recovers spatial beating and rotation, and identifies a
coarse-mesh equilibrium absent from the continuum model. Its stronger-field
scout fails the sampled clearance screen. The planar records on this page
remain unchanged.

The 2026 article, its 2024 planar predecessor and targeted author/model/chaos
queries were checked. This is a bounded source audit, not a comprehensive
citation or equivalence search. Absence of a matching search result does not
establish historical novelty. Chaos in other active-filament models is already
known; the earlier [follower-force screen](geometry-search-20261031.md) also
encountered self-intersections.

The unresolved target is the full three-dimensional morphoelastic system,
including twist equilibrium and transverse stability. A spatial solver must
recover these planar controls and the source's spatial periodic examples
before its complicated trajectories can be used as evidence. The present
packet closes neither that numerical task nor the search for a new attractor.

[Research index](../README.md#research-index) ·
[Candidate queue](candidate-queue.md) · [Reproduction](reproduction.md)

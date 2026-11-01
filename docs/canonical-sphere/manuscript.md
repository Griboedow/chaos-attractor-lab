# Folded ribbons and spatial networks in canonical posterior-reflection maps

**Numerical research article · 2026-11-01.** Asymptotic chaotic attraction
and historical uniqueness remain unresolved.

## Abstract

We simplify a Gaussian posterior-reflection iteration to regular tetrahedral
and cubic sites with short decimal parameters. The ideal posterior is exactly
a three-spin mean: independent spins for a cube and one parity constraint
for a tetrahedron. A complete 656-case exploration produces a folded ribbon
and a spatial network of curved sheets. Both selected regimes retain positive
leading growth over five million iterations, 33 fresh starts, and independent
70- and 100-digit calculations. The cubic candidate also has two positive
finite-time QR rates. Periodic numerical captures, including unstable fixed
points with positive tangent growth, are excluded from promotion. The
contribution is a reproducible numerical parameter study and explicit
equivalence accounting. Neither polyhedral appearance, sigmoid feedback,
nor hyperchaos is a new general mechanism.

## 1. Visual candidates

![Canonical posterior maps, Cases 163 and 653](../assets/portraits/canonical-sphere.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Case 163 forms an elongated ribbon with several curled edges; Case 653
forms a spatial network of folded sheets. Each panel contains two million
consecutive states of a native three-coordinate discrete map. These are
density images, not interpolated curves or rendered surfaces. The
[24-second rotation](../assets/portraits/canonical-sphere.webp) keeps the
terminal 250,000 states fixed and rotates only the camera. Centering and
orthogonal rotations preserve geometry. Units are equal within each panel;
panel scales and density saturation are independent.

**Case numbers are local to this canonical-site campaign.** They do not
refer to the identically numbered cases of earlier studies. The
[preceding sphere article](../soft-sphere/manuscript.md) uses rotated,
separately rounded site arrays and a different temperature sweep.

## 2. Model and exact ideal reductions

For sites $`b_\sigma=c+a\sigma`$, define

```math
p_\sigma(x)=\frac{\exp(-\|x-b_\sigma\|^2/\tau)}
 {\sum_\eta\exp(-\|x-b_\eta\|^2/\tau)},\quad
m(x)=\sum_\sigma p_\sigma(x)b_\sigma,\quad
F(x)=\beta x+\beta\frac{2m(x)-x}{\|2m(x)-x\|}+(1-2\beta)m(x).
```

The sphere is centered at the origin and has radius one. For cubic sites,
$`\sigma\in\{-1,1\}^3`$. For tetrahedral sites retain the four signs with
$`\sigma_1\sigma_2\sigma_3=1`$. The map is smooth wherever
$`2m(x)-x\ne0`$. The derivative, reflection interpretation and elementary
boundedness inequality are given in the
[preceding article, Section 2](../soft-sphere/manuscript.md#2-map-derivative-and-elementary-bound).
Boundedness does not prove avoidance of the singular set or chaotic attraction.

Put $`h=2a(x-c)/\tau`$ and $`t_i=\tanh h_i`$. Factorizing the cube's
partition function and restricting it to even parity for the tetrahedron gives

```math
m_i(x)=c_i+a t_i\quad\text{(cube)},\qquad
m_i(x)=c_i+a\frac{t_i+t_jt_k}{1+t_1t_2t_3}
\quad\text{(tetrahedron)},\quad \{i,j,k\}=\{1,2,3\}.
```

Thus the four-site model is also an exact posterior for a three-bit
even-parity code, with signs $`\sigma_i=(-1)^{q_i}`$. This is an algebraic
interpretation of the posterior, not an identity between the full iteration
and belief propagation. The production implementation retains shifted
exponentials: direct evaluation of the tetrahedral quotient can produce
`0/0` when floating-point tangents saturate.

Orthogonal coordinate changes preserve the sphere and remove the initial
site rotation. Translation alone does not. The identities above hold for
ideal regular sites. Executed inputs are separately rounded binary64 arrays
computed from the displayed parameters; their exact hexadecimal values are
retained. Numerical reconstruction of an earlier rotated array is not an
exact conjugacy certificate for that array. The
[shared equivalence review](../soft-constellation-priority.md#canonical-spin-and-neural-formulations)
gives the derivation and related neural formulations.

| Case | Site set | beta | tau | a | Center c |
| --- | --- | --- | --- | --- | --- |
| 7 | Tetrahedron | 0.93 | 0.059 | 0.8 | (-0.1, -0.2, 0.2) |
| 148 | Tetrahedron | 0.93 | 0.06 | 0.8 | (-0.1, -0.2, 0.2) |
| **163** | **Tetrahedron** | **0.93** | **0.1** | **0.8** | **(-0.1, -0.2, 0.2)** |
| 361 | Tetrahedron | 0.99 | 0.2 | 0.8 | (0.1, 0.1, 0.1) |
| 407 | Cube | 0.9 | 0.02 | 0.25 | (0, 0, 0) |
| 565 | Cube | 0.96 | 0.04 | 0.3 | (-0.5, -0.5, -0.5) |
| 638 | Cube | 0.99 | 0.04 | 0.3 | (-0.6, -0.4, -0.2) |
| **653** | **Cube** | **0.99** | **0.08** | **0.35** | **(-0.5, -0.5, -0.5)** |

## 3. Complete selection history

The screen contains eight reconstructions or joint decimal simplifications
of the earlier sphere Cases 155 and 251, followed by complete Cartesian
grids of 360 tetrahedral and 288 cubic inputs. The
[protocol](../../results/canonical-sphere/runs/canonical-sphere-scout01/protocol.json)
lists every parameter and selection reason. Every run starts at
`(0.317,-0.211,0.419)`. All 656 complete a 5,000-step discard and 15,000
observations. The 262 cases whose leading rate exceeds 0.005 in every one
of ten windows receive fresh 30,000-step discards and 500,000 observations.

Of those follow-ups, 254 retain positive leading windows and 18 have
detected periods from 1 to 512. Ten satisfy both conditions. Excluding every
period flag leaves 244 positive, nonrepeating numerical candidates. These
are finite tests with a scale-adjusted `1e-8` terminal return tolerance,
not claims of asymptotic aperiodicity or 244 distinct attractors.

All 22 pages of the [complete preview](../assets/experimental/canonical-sphere-first-look/render.json)
were visually inspected. Cases 7, 163, 361, 407, 565, 638 and 653 were
selected for fresh runs discarding 100,000 iterations and observing five
million in twenty windows. Case 148 was subsequently added to test the
simplest neighboring temperature, 0.06, rather than 0.059. Its selection
is explicit and retained separately. All eight outcomes are reported.

| Case | Mean finite-time QR rates, per iteration | Smallest leading window | Terminal period flag |
| --- | --- | --- | --- |
| 7 | 0.70302, -0.01148, -0.84507 | 0.70033 | None through 512 |
| 148 | 0.66329, 0.00750, -0.84556 | 0.65990 | None through 512 |
| **163** | **0.20665, -0.04250, -0.46071** | **0.20497** | **None through 512** |
| 361 | 0.68171, -0.32731, -1.04273 | 0.67874 | None through 512 |
| 407 | 0.28220, 0.12535, -0.04682 | -0.05151 | 6 |
| 565 | 0.19000, -0.08930, -0.98743 | 0.18865 | None through 512 |
| 638 | 0.16237, -0.33747, -1.93404 | 0.16111 | None through 512 |
| **653** | **0.16191, 0.01395, -0.62918** | **0.15973** | **None through 512** |

Case 407 is rejected despite its initially elaborate spherical portrait:
it enters a detected six-cycle. Its positive run average includes a long
transient; later windows are negative. Case 638 has a smallest covariance
fraction of only `3.58e-5` and is nearly planar. For Cases 163 and 653 the
three terminal covariance fractions are respectively
`(0.01764,0.03086,0.95150)` and `(0.10210,0.41866,0.47925)`.
The ribbon is thin, whereas the cubic network has substantial extent in
all three principal directions.

## 4. Precision, fresh starts and false positives

For Cases 7, 163, 361 and 653, the original start and 32 scrambled Sobol
starts in `[-4,4]^3` each receive a 100,000-step discard and 500,000
observations. All 33 starts for each case retain positive leading growth
in every window. For the two displayed cases the minimum windows over
all starts are 0.20194 and 0.15482. Finite basin sampling does not establish
global attraction or a unique invariant measure.

None of these 132 cohort trajectories has a terminal period flag through
512. Iterating the posterior mean alone reaches a `1e-12` stopping residual
from all cohort starts, in at most 3, 4, 4 and 10 steps for Cases 7, 163,
361 and 653. This control is consistent with known mean-shift convergence;
it does not transfer a descent property to the reflected composition.

Independent scalar calculations use posterior inner-product logits and
nested reflections, without the production Jacobian. Exact binary64
inputs are lifted to 70 and 100 decimal digits. After 10,000 discarded
steps, renormalized finite separations are observed for 50,000 steps;
separations are `1e-30` and `1e-40`, respectively. Case 163 has minimum
leading windows 0.20314 and 0.19974; Case 653 has 0.15207 and 0.15784.

Case 653 has two positive QR rates in all twenty 250,000-step windows of
the long run; the smallest second rate is 0.01234. A separate calculation
propagates three finite separations and orthogonalizes them in high precision.
At 70 digits it gives `(0.16088,0.01401,-0.62008)` over 100,000 observations,
with smallest second window 0.00511. This supports a numerical hyperchaos
candidate; it is not an interval proof of positive asymptotic exponents.
At 100 digits the mean spectrum is `(0.16179,0.01372,-0.62585)`, but the
smallest second 5,000-step window is negative, `-0.00963`. Agreement of
the positive long averages therefore does not imply uniformly positive
short-time expansion.

Full QR calculations from the same 33 starts also retain two positive
rates in every 25,000-step window for Case 653; the smallest second window
is 0.00419. Case 148 is less uniform: all 33 full-run second means are
positive, but only four starts have positive second rates in every shorter
window. Its 70-digit second mean is 0.00729 with negative short windows.
The 100-digit second mean is 0.00427, also with negative short windows.
The full window records are retained rather than replacing this variation
by a single favorable average.

The rejected short-run Cases 35 and 300 illustrate a different failure mode.
Their stored terminal trajectories are exactly constant in binary64, yet
their computed largest multipliers are about 19.73 and 128. Independent
80-digit trajectories starting `1e-40` apart in the computed unstable
direction separate by more than 0.01 after 30 and 19 iterations. A positive
tangent rate at a numerically captured unstable fixed point is not evidence
of a chaotic attractor. The [inspection record](../../results/canonical-sphere/runs/canonical-sphere-inspection01/report.json)
preserves these controls; neither a certified root nor a certified basin is
inferred from them.

## 5. What is established, and what was already known

The reduction exposes established descriptions of the same averaging step:
Gaussian mean shift, finite-prior MMSE denoising, softmax attention, and
spin posteriors. Leaky sigmoid chaos appears in
[Chen and Aihara (1997)](https://arxiv.org/html/chao-dyn/9701020v1);
radial neural normalization appears in
[Verzelli, Alippi and Livi (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6761167/).
Neither reviewed recurrence is literally the displayed posterior-reflection
map. That formula distinction does not exclude a nonlinear conjugacy.

The parity interpretation likewise does not establish a first chaotic
decoder: [Sibel, Reynal and Declercq](https://perso.etis-lab.fr/declercq/PDF/ConferencePapers/Sibel_2012_CMSIM.pdf)
report chaos in iterative belief propagation. Their state consists of graph
messages, whereas this study reflects a three-coordinate exact posterior.
The [shared review](../soft-constellation-priority.md) records which equations
and assumptions were actually read.

Polyhedral chaotic maps predate this study by decades. Reiter's
[n-cube maps (1996)](https://emis.de/ft/51972),
[tetrahedral study (1997; author abstract)](https://webbox.lafayette.edu/~reiterc/abstracts/mv_1997c.html),
and [dodecahedral maps (1999)](https://webbox.lafayette.edu/~reiterc/dodec/dodec02_pp.pdf)
are direct comparison targets. The last paper's Figures 1–2 already display
spatial threads, loops and diffuse interior structure. Its polynomial
group-averaged construction differs from our posterior-reflection formula,
but broad visual motifs and symmetry cannot establish uniqueness.

Two positive exponents in a three-dimensional discrete map are also known.
For example, [Toker, Sommer and D'Esposito (2020), equation (12)](https://www.nature.com/articles/s42003-019-0715-9)
uses the generalized Hénon recurrence

```math
s_{n+1}=1.76-s_{n-1}^2-0.1s_{n-2}
```

as a hyperchaotic reference. Thus the sign pattern of Case 653's estimated
spectrum is a useful dynamical observation, not a first-discovery claim.

A separate [finite cloud comparison](../../results/canonical-sphere/runs/canonical-sphere-catalog01/report.json)
uses 108 frozen reference regimes. Under similarity normalization, the
first fingerprint-ranked references are `WindmiReduced` for Case 163 and
`SprottM` for Case 653; the fitted local ICP residuals are 0.19183 and
0.24189, against time-half self-residuals 0.02489 and 0.04557. After affine
whitening, the first labels become `BelousovZhabotinsky` and `SprottQ`.
These are archived software labels, not newly verified reproductions of
every cited primary model. Only three fingerprint-ranked references per
representation receive 48-start local alignment. The comparison neither
finds a global geometric optimum nor tests nonlinear conjugacy, other
parameter regimes, or every known attractor.

The present result is a short-parameter, reproducible pair of numerical
geometries, together with exact ideal posterior reductions and explicit
negative controls. The full iteration is not proved equivalent to the
reviewed examples; neither is it proved inequivalent to them. Certified
attraction, invariant structure, nonlinear equivalence and historical
priority remain open. The certified theorems for earlier circle maps do
not transfer to these sphere candidates.

## 6. Reproduction

The [evidence packet](../../results/canonical-sphere/README.md) retains the
complete search, source snapshots, exact inputs, numerical controls and
commands. Large unthinned arrays remain in a hash-addressed external cache.
Replay checks account for every screen output and follow-up, repeat the
long observations and dense rendering, and compare local QR arithmetic
with LAPACK and determinant accumulation. These checks establish numerical
reproducibility in the recorded environment, not a mathematical chaos proof.

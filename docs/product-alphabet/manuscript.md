# Folded wings in product-alphabet posterior-reflection maps

**Computational research article · 2026-11-01.** Positive entropy is certified
on a compact invariant subsystem. Chaotic attraction and historical uniqueness
are not certified.

## Abstract

An exact Gaussian posterior for three independently quantized coordinates,
combined with unit-sphere reflection feedback, produces a spatial network of
four folded wings. A declared 1,753-case survey covers binary, ternary and
four-level Cartesian alphabets. All 175 initial nonperiodic flags receive
longer follow-up; twelve visually selected regimes receive five-million-step
observations. The principal binary candidate retains positive leading growth
at all 33 tested starts and in independent 70- and 100-digit calculations.
A separate interval certificate proves a compact invariant subsystem with
topological entropy above 0.0254 per original iteration for the exact binary
parameters. A second, full-eight-site automatic-derivative implementation
verifies every covering relation. This subsystem is not proved attracting
or equal to the sampled cloud. A four-level companion produces layered arches
but coexists with numerical fixed points. Complete selection records preserve periodic captures and
precision-sensitive planar behavior. The contribution is a reproducible
geometric study and an explicit positive-entropy certificate. Posterior means,
spin responses, reflection algorithms and feedback-induced chaos are known;
neither a new general operator nor historical priority is claimed.

## 1. Geometry and evidence level

![Product-alphabet posterior maps, Cases 988 and 489](../assets/portraits/product-alphabet.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Case 988 is the principal candidate: four curved, folded wings occupy a
substantial volume in all three coordinate directions. Case 489 is a
more anisotropic companion with layered arches and curled edges. These
descriptions refer to numerical point clouds, not certified surfaces or
topological classifications. Each panel includes two million consecutive
states after 100,000 discarded iterations. No symmetry replication,
interpolation, smoothing, or nonlinear visual transformation is applied.
Centering and orthogonal cameras preserve native geometry; spatial units
are equal within each panel. Panel scales and density saturation differ.

The [24-second rotation](../assets/portraits/product-alphabet.webp) uses
the last 250,000 consecutive states of the same rendering inputs, with
fixed states and scale throughout. The apparent curves are dense sets of
discrete iterates. No period at most 512 is detected in the recorded
selected tails; this finite test does not prove aperiodicity.

**Case identifiers are local to this product-alphabet campaign.** The
[canonical-site study](../canonical-sphere/manuscript.md) has separate case
numbers. In particular, this Case 988 is not a renaming of its cubic Case 653.
The full [evidence packet](../../results/product-alphabet/README.md) retains
exact inputs, source snapshots, unsuccessful selections and numerical checks.

## 2. Model, factorization and established interpretations

For coordinate alphabets $`A_i=\{\ell_{i,0},\ldots,\ell_{i,q-1}\}`$,
take the uniform prior on $`A_1\times A_2\times A_3`$. Independent Gaussian
observation noise with variance $`\tau/2`$ gives

```math
w_{i,k}(x_i)=\frac{\exp(-(x_i-\ell_{i,k})^2/\tau)}
 {\sum_{j=0}^{q-1}\exp(-(x_i-\ell_{i,j})^2/\tau)},
\qquad m_i(x_i)=\sum_{k=0}^{q-1}w_{i,k}(x_i)\ell_{i,k}.
```

The full posterior factors exactly into these three scalar distributions.
The study evaluates $`q=2,3,4`$, corresponding to 8, 27 and 64 spatial sites.
The deterministic recurrence is

```math
v(x)=2m(x)-x,\qquad r(x)=\|v(x)\|,\qquad u(x)=v(x)/r(x),
\qquad F(x)=\beta x+\beta u(x)+(1-2\beta)m(x).
```

It is smooth wherever $`r(x)>0`$. If $`R_m=2m-I`$ and
$`R_S(v)=2v/\|v\|-v`$, then

```math
F=\frac{\beta}{2}(R_S R_m+I)+(1-\beta)m.
```

This is the established relaxed-reflection construction with a posterior
mean in place of one metric projector. The signal interpretation concerns
this deliberately specified iteration, not an observed communication
device or a standard decoder implementation. The sphere introduces a
known coupling through normalization. The [shared equivalence review](../soft-constellation-priority.md)
covers mean shift, MMSE estimation, proximal interpretations, plug-and-play
splitting, sigmoid networks and chaotic iterative decoding.

Write $`D=\mathrm{diag}(d_1,d_2,d_3)`$. The exact derivative is

```math
d_i=\frac{2}{\tau}\sum_k w_{i,k}(\ell_{i,k}-m_i)^2,
\qquad
DF=\beta I+(1-2\beta)D+
\frac{\beta}{r}(I-uu^{\mathsf T})(2D-I).
```

For an ideal symmetric binary alphabet $`\ell_{i,\pm}=c_i\pm a_i`$,
the familiar independent-spin response is

```math
h_i=\frac{2a_i(x_i-c_i)}{\tau},\qquad
m_i=c_i+a_i\tanh h_i.
```

For the ideal ternary alphabet $`c_i+a_i\{-1,0,1\}`$, it becomes

```math
m_i=c_i+a_i\frac{2\sinh h_i}{\exp(a_i^2/\tau)+2\cosh h_i}.
```

The latter is the standard spin-1/Blume–Capel single-site response, after
parameter identification. These identities prevent presenting the scalar
nonlinearity as a new invention. Executed alphabets are separately rounded
binary64 arrays built from $`c_i+a_i\,\mathrm{linspace}(-1,1,q)`$.
Their exact entries, rather than ideal symmetry identities or rounded
table entries, define the recorded experiments.

| Case | Levels per coordinate | beta | tau |
| --- | --- | --- | --- |
| **988** | **2** | **0.9268088431125507** | **0.020935697765931374** |
| 489 | 4 | 0.9699380761357025 | 0.035825792302107334 |

| Case | Halfwidths a | Center c |
| --- | --- | --- |
| 988 | (0.9987935485318303, 0.9138440923765302, 1.1366534527391194) | (-0.13522116094827652, 0.2182737350463868, -0.07052738964557648) |
| 489 | (1.0288301357999443, 0.17564048673957586, 0.18828875999897718) | (0.104424723982811, -0.5549065947532654, -0.0760837271809578) |

Exact hexadecimal parameters and alphabet entries are recorded in the
verification report. Translation or unequal coordinate scaling does not
preserve the unit-sphere normalization, so simply centering or stretching
the alphabet does not produce an immediate conjugacy to the earlier cubic
map. This observation excludes only that elementary identification.

## 3. Complete search and selection

The initial state is $`(0.317,-0.211,0.419)`$. The declared screen contains
one exact input control from canonical Case 653, a complete 216-case
centered isotropic grid, and 512 scrambled Sobol draws each evaluated at
all three alphabet orders. The Sobol seed is `2026110144`. The eight
coordinates determine beta, temperature, three halfwidths and three shifts;
the exact ranges and construction appear in the executed source.

After 5,000 discarded and 15,000 observed iterations, 177 cases have a
leading rate above 0.005 in all ten windows. A period test through lag 64
rejects two. All remaining 175 receive 30,000 discarded and 500,000
observed iterations with full three-direction QR calculations. Of these,
165 retain positive leading growth in every window without a detected
period: 46 binary, 60 ternary and 59 four-level cases. Nine others have
detected periods (four fixed, three period two, one period four and one
period 43); one additional case has a nonpositive leading window.

All fifteen preview pages were inspected before selecting twelve cases:
203, 287, 397, 489, 607, 893, 933, 988, 1351, 1450, 1660 and 1709.
Each receives a fresh 100,000-step transient, five million observed steps,
twenty full-spectrum windows, and a separate two-million-state portrait.
All twelve retain positive leading growth. Eight have an observed smallest
covariance fraction below $`10^{-5}`$; numerical growth alone does not make
them convincing spatial candidates. Five cases receive the start and
precision controls below. The two displayed cases were selected for
geometry and those controls, not through a preregistered visual score.

## 4. Numerical validation and rejected interpretations

| Diagnostic | Case 988 | Case 489 |
| --- | --- | --- |
| Five-million-step QR rates per iteration | (0.251100, -0.280399, -0.440762) | (0.337780, -0.406624, -0.787602) |
| Smallest leading window | 0.248311 | 0.335313 |
| Covariance fractions, ascending | (0.237213, 0.327906, 0.434880) | (0.020710, 0.086140, 0.893149) |
| Positive, nonperiodic outcomes among 33 starts | 33 | 13 |
| Numerical fixed-point outcomes among 33 starts | 0 | 20 |
| Smallest 70-digit leading window | 0.235852 | 0.319821 |
| Smallest 100-digit leading window | 0.229787 | 0.322401 |

The start cohort combines the original initial condition and 32 scrambled
Sobol points in $`[-4,4]^3`$, seed `2026110146`. Each gets 100,000 discarded
and 500,000 observed iterations, twenty full-spectrum windows and a period
test through lag 512. Case 988 has a minimum leading window of 0.232913
over the entire cohort. The finite sample does not establish basin volume,
a unique attractor, or a global basin. Case 489's coexistence is a material
limitation. Other tested cases are retained: Case 933 gives seven positive
outcomes, twenty fixed points and six periodic captures; Case 1450 gives
eleven positive outcomes and twenty-two fixed points.

The arbitrary-precision controls lift the exact binary64 inputs into
70- and 100-digit arithmetic. They use scalar squared-distance Gaussian
weights, nested reflections and finite separation, without the production
Jacobian. Each discards 10,000 iterations and observes 50,000 in ten windows;
separation sizes are $`10^{-30}`$ and $`10^{-45}`$, respectively. This is an
independent numerical expression and precision control, not an interval
proof or an independent derivation of posterior factorization.

A separate 70-digit check sums all eight spatial posterior weights for
Case 988, without scalar factorization. Its mean leading rate is 0.245820
and its smallest window is 0.223730 over the same observation length.
This additional implementation retains the nested-reflection expression
and finite-separation method, so those components are shared.

Parameter neighborhoods contain regular windows. For each selected case,
64 Sobol directions are evaluated at two declared perturbation radii.
At radius 0.01, 45 of 64 runs retain positive leading growth without a
detected period for each case. At radius 0.05, the corresponding counts
are 39 for Case 988 and 26 for Case 489; all other runs in these cohorts
have detected periods. Beta changes by the radius times its distance from
one, temperature and halfwidths change multiplicatively in logarithmic
coordinates, and centers change in halfwidth units. Each run uses the
original initial state, a 100,000-step transient and 500,000 observed
iterations. These results support nearby recurrence but do not establish
an open region of robust chaos or preservation of the four-wing geometry.

**Case 203 exposes a precision-sensitive capture.** Its binary64 dense
trajectory lies on the exact floating plane $`x=z`$ and has leading growth
about 0.809. The high-precision calculations instead retain three comparable
covariance fractions and leading growth about 0.56 over their recorded
windows. Its attractive-looking projected cloud therefore cannot be
reported as evidence for the same spatial regime across precisions. No
claim about the infinite-time fate at arbitrary precision follows.

All screen arrays, all 175 follow-ups, all twelve long observations and all
twelve dense trajectories and histograms are replayed from the stored
inputs. Independent full-site posterior/covariance evaluations and nested
finite differences check the factorized map and derivative at 192 points.
QR steps are checked against LAPACK with the same incoming frame, and
accumulated volume changes against log determinants. Independently propagated
subordinate QR frames differ in two nearly planar cases; the first verifier
rejected this discrepancy and is preserved. The corrected verifier records
both rates and tests the local identity separately. Both selected candidates
also pass the accumulated three-rate agreement check. Numerical replay uses
the production recurrence and does not independently certify its dynamics.

## 5. Approach from multiple starts

![Thirty-three initial states approaching the Case 988 cloud](../assets/diagnostics/product-alphabet-approach.png)

The left panel shows the broad initial distribution; the middle panel
enlarges the central region. Colored dots are actual iterates 0 through
100, with color indicating iteration. The right panel shows the median
and 10–90 percentile band of the distance to a separate 100,000-state
reference tail, normalized by that tail's bounding-box diagonal. The
distance is a finite-cloud diagnostic, not a certified distance to an
invariant set. It need not approach zero at the displayed resolution.
The central zoom clips the distant initial points explicitly shown on the
left. No continuous orbit is inferred between dots.

For comparison, iterating the posterior mean alone from the same 33 starts
reaches the numerical stopping threshold within three steps in Case 988
and 82 steps in Case 489. Thus these sampled irregular regimes depend on
the reflection feedback. This control does not claim that posterior-mean
iterations are globally convergent for every finite prior.

## 6. A certified positive-entropy subsystem

The [separate symbolic packet](../../results/product-symbolic/README.md)
uses exactly the binary64 parameters and eight sites of Case 988. Put
$`G=F^4`$. Twenty-three affine windows
$`N_i=c_i+S_i[-1,1]^3`$ have one exit coordinate and two entry coordinates.
The centers, basis factors and radii are exact real lifts of stored binary64
values. The matrices $`S_i`$ and their inverses are formed over the rationals;
their complete entries are retained in
[input.json](../../results/product-symbolic/runs/product-coverings988-certificate04/input.json).

The certified graph is

```math
N_0\longrightarrow N_0,\qquad
N_0\longrightarrow N_1\longrightarrow\cdots\longrightarrow N_{22}
\longrightarrow N_0,
\qquad\text{each arrow represents }G=F^4.
```

For each arrow, the full source volume maps either strictly outside the
target exit interval, or strictly inside both target entry intervals.
The two exit faces map strictly to opposite exterior sides. These are
sufficient conditions for Theorem 16 of
[Zgliczynski and Gidea (2004)](https://ww2.ii.uj.edu.pl/~zgliczyn/papers/wts/cov.pdf).
Theorem 9 of that work realizes every closed chain by a periodic point.
These established covering theorems supply the topological step; the new
calculation supplies their hypotheses for the exact Case 988 map.

**Two interval implementations.** The primary calculation uses 80-digit
outward intervals, the exact binary posterior factorization and its analytic
Jacobian. It composes four original steps and applies mean-value enclosures
in exact rational window coordinates. The second calculation uses 110-digit
intervals, all eight Gaussian squared-distance weights, nested spherical
reflections and forward automatic derivatives. Its weighted sum is centered
at one site by an exact algebraic identity to reduce interval overestimation.
It revisits every primary leaf in reverse order and subdivides further when
needed. Both implementations share mpmath interval arithmetic, the supplied
windows and the mean-value enclosure principle.

| Coverings | Primary leaves | Full-site AD leaves |
| --- | ---: | ---: |
| 0 to 0 | 16 | 16 |
| 0 to 1 | 11 | 11 |
| Remaining 22 arrows | 151 | 183 |
| Total, including volumes and both exit faces | 178 | 210 |

Exact-rational audits check complete binary subdivision trees, all strict
inequalities using stored binary interval endpoints, matrix inverses,
source hashes and $`N_0\cap N_1=\varnothing`$. Minimum normalized margins
exceed 0.00629 and 0.01915, respectively. The sphere-normalization radius
exceeds 2.3799 throughout all four-step source enclosures in both checks.
These are bounds on the validated windows and their intermediate images,
not on the entire numerical cloud or a parameter neighborhood.

![Locations of covering-window centers and their two return routes](../assets/diagnostics/product-alphabet-coverings.png)

**Figure 3.** The background is the same two-million-state numerical cloud
seen through an orthogonal, equal-unit camera. Orange markers locate window
centers; the blue diamond marks N_0. Symbols are enlarged, and many centers
overlap at this scale. The right panel is schematic: 23 self-loops and one
23-arrow excursion each use 92 original iterations. Window centers are not
asserted to constitute an exact orbit. Attraction is not inferred from
their plotted proximity to the cloud.

**Theorem.** For the exact Case 988 instance, there is a nonempty compact set
$`C`$ in the smooth domain with $`F(C)=C`$ and

```math
h_{\mathrm{top}}(F|_C)\ge\frac{\log\rho}{4}>0.0254,
\qquad \rho>1,\qquad \rho^{23}-\rho^{22}-1=0.
```

**Proof.** A return to N_0 may use either one self-loop or the 23-arrow
excursion. Let $`a_n`$ count words in these two returns with total length n
in G-steps. Then $`a_0=1`$, $`a_n=0`$ for negative n, and
$`a_n=a_{n-1}+a_{n-23}`$. Each word gives a closed covering chain and hence
a point fixed by $`G^n=F^{4n}`$ with its prescribed visits.

Let $`\delta=\mathrm{dist}(N_0,N_1)>0`$. Two distinct equal-length words have
a first different return choice after a common prefix. One G-step later
their realizing points lie in N_0 and N_1. Therefore these $`a_n`$ points
are $`(4n,\delta/2)`$-separated for F. This argument does not require the
other windows to be pairwise disjoint.

All these periodic orbits lie in the compact set

```math
X=\bigcup_{i=0}^{22}\ \bigcup_{j=0}^{3}F^j(N_i).
```

The validated normalization bound puts X in the smooth domain. Let C be
the closure of the union of the selected periodic orbits. Continuity gives
$`F(C)\subset C`$. Compactness of X supplies a convergent subsequence of
periodic predecessors for each limit point, giving $`C\subset F(C)`$.
Thus $`F(C)=C`$. The positive root rho above obeys
$`\rho^n=\rho^{n-1}+\rho^{n-23}`$. Upper and lower positive multiples of
$`\rho^n`$ bound $`a_n`$, first for n from 0 to 22 and then by induction.
The separated-set definition of entropy gives the stated bound. The audit
checks $`(1.10695)^{22}(1.10695-1)<1`$ exactly and bounds its logarithm
by a positive rational series, establishing the displayed decimal lower
bound without a floating-point root assumption. ∎

As a simpler consequence, choosing 23 self-loops or one excursion in each
block gives arbitrary binary itineraries sampled every 92 original steps
and the weaker bound $`\log(2)/92>0.0075342`$. Compactness also realizes
infinite prescribed routes. Unique symbolic coding, a transverse homoclinic
orbit, uniform hyperbolicity and attraction of C are not asserted.

**Search and unsuccessful certificates.** The window search starts from
512 numerical cycle fits at periods 1–16, retaining 70 numerical cycles.
A period-four saddle proposal guides linear unstable-line searches and
multiple shooting with linear endpoint constraints. These guides are not
used as invariant-manifold certificates. Early window programs failed;
the default absolute linear-program feasibility tolerance also exceeded
some proposed widths. Rescaling and explicit constraint checks corrected
that numerical issue before the interval validation. Two grid-trapping
proposals did not close. A wider, uncentered full-site interval expression
was stopped incomplete after the algebraically centered implementation
passed. All these outcomes and executed sources are inventoried separately
from the frozen numerical survey. None of the unsuccessful trapping runs
is used to claim that C attracts a neighborhood.

## 7. Prior art and geometric comparison

The literature audit explicitly includes equivalent terminology. Its
[receipt](../../results/product-alphabet/literature-review.json) distinguishes
reviewed equations, abstract-only evidence and uninspected figures.

| Earlier formulation | Relationship and limit |
| --- | --- |
| MMSE/mean shift, proximal denoisers, relaxed reflections | Established ingredients of the full construction; reviewed in the shared priority note. |
| Independent-spin tanh response | Exact ideal binary posterior identity. Calling it a new quantizer or neural activation would be incorrect. |
| Spin-1/Blume–Capel response | Exact ideal ternary posterior identity. The scalar response is established. |
| Chaotic spin-1 recursion | Tomé and Salinas (1989) already report strange-attractor evidence in a two-dimensional competing-interaction recursion. |
| Spin models with feedback | Sinelschikov and colleagues (2023) study oscillations and chaos in feedback lattice models, including a four-dimensional mean-field Potts flow. |
| Labyrinth and polyhedral attractors | Spatial loops, folded threads and polyhedral organization have direct precedents; their appearance alone is not novel. |
| Iterative phase retrieval | RAAR and RRR trajectories can retain stationary fluctuations; this algorithmic behavior predates the present examples. |

[Ülgen and Sünel (2013)](https://dergipark.org.tr/en/download/article-file/192677)
give the Blume–Capel mean-state recurrence in Section 2.1. The reviewed
scalar formula matches the ideal ternary response after parameter
identification, but does not contain the present three-coordinate radial
feedback. [Tomé and Salinas (1989)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.39.2206)
report chaos in a spin-1 recursion; only their publisher abstract was
accessible in this audit, so no full equation or figure comparison is claimed.

[Sinelschikov, Poggialini, Abbate and De Martino (2023)](https://arxiv.org/abs/2306.01823)
derive mean-field feedback dynamics for Blume–Capel and Potts systems.
Sections 2.3–2.4, including equations (14)–(34) in the reviewed preprint,
establish a particularly relevant precedent for softmax/spin feedback and
chaos. Their reduced Potts model is a four-dimensional continuous flow,
whereas the present recurrence has three discrete coordinates and explicit
sphere normalization. This difference in displayed equations is not a
proof against nonlinear equivalence or equivalence to a return map.
Figure 5 was also visually inspected: its chaotic panels show Poincare
sections in two coordinates, so their curled islands are not directly
comparable to a full three-state cloud. The inspected figure and surrounding
text do not supply an equation-level identification with Case 988.

The finite geometric comparison uses 108 frozen reference regimes.
Fingerprint ranking is followed by local point-cloud alignment for the
three highest-ranked references, using 48 initial orthogonal alignments
and twenty local iterations on 1,500 points. Similarity normalization and
affine whitening are tested separately. For Case 988, the closest
similarity fingerprint is the frozen `ThomasLabyrinth` trajectory; its
local alignment residual is 0.24454, compared with 0.04702 between two
time halves of the candidate. After affine whitening, the leading
fingerprint is `Coullet`, with local residual 0.22509 versus a self-comparison
of 0.05269. These are dataset/code labels, not independently revalidated
historical attributions for every reference formula.

[Sprott and Chlouverakis, *Labyrinth Chaos* (2007)](https://sprott.physics.wisc.edu/pubs/paper302.pdf)
provide a primary precedent for the nearest similarity family. Its cyclic
sine flow, coexistence discussion and Figures 4–9 were reviewed. Figure 8
combines six coexisting strange attractors, whereas each of our selected
portraits comes from one orbit. This specific distinction does not rule out
similar single-attractor geometry elsewhere in that family. Case 489's closest similarity fingerprint
is `Chua` (local residual 0.16861 versus self-residual 0.04135). All 108
rankings and both normalization controls are retained, including references
whose local fit is better despite a worse fingerprint rank. No global
alignment optimum, topological invariant, nonlinear conjugacy test or
historical exclusion follows from these finite distances.

[Barnett, Millane and Kingston (2024)](https://journals.iucr.org/d/issues/2024/11/00/nz5017/index.html)
compare DM, RRR and RAAR in crystallographic phase retrieval. Their Section
2.1.4, equation (6), is the same RAAR algebra with metric projections;
Sections 6–7 discuss persistent fluctuations around a stationary mean.
Those observations prevent claiming the first irregular reflected
reconstruction trajectory. The reviewed constraints are Fourier amplitudes,
solvent flatness and histogram matching, rather than this three-coordinate
product posterior. The text and Figure 7 caption were reviewed; the image
download was inaccessible, so no visual comparison to that figure is claimed.
The supplementary [review receipt](../../results/product-symbolic/literature-review.json)
records this scope and the covering-theorem source.

## 8. Contribution and open questions

The supported result is a fully specified four-wing numerical regime of
an eight-site product posterior with sphere feedback, a multilevel
companion with documented coexistence, and a complete search and control
record. A separate certificate proves positive entropy on a compact invariant
subsystem for the exact principal parameters. The unexpected setting is a quantized-signal posterior inside an
established numerical reflection construction. It is not a new physical
phenomenon established by measurement.

Neither a strange-attractor theorem nor a claim that the geometry has
never appeared before is justified yet. Positive finite-time growth,
multiple starts, precision controls and finite-library differences narrow
the candidates; those numerical tests do not prove asymptotic chaos,
transitivity, a fractal dimension, or historical uniqueness. Section 6 supplies
an independent positive-entropy theorem, but leaves attraction and equality
with the depicted cloud open. A useful next theoretical target is a verified
trapping region containing the certified subsystem, followed by structural
comparison with known noninvertible maps. The
independent certificates for earlier cases do not transfer to this parameter
set merely because the formulas belong to the same family.

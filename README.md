# Chaos Attractor Lab

Computational studies of unusual attracting dynamics in continuous systems and numerical algorithms. The repository separates **certified dynamical properties**, **numerical chaos candidates**, and **unresolved or negative results**. Historical novelty is not established by visual appearance or by an unsuccessful literature search.

[Portrait gallery](docs/attractor-gallery.md) · [Research index](#research-index) · [Manuscripts](papers/README.md) · [Reproduction](docs/reproduction.md) · [Data](results/README.md) · [Storage](docs/data-storage.md)

The portraits rotate once every 24 seconds at 25 frames per second. Each view uses fixed trajectory data and a constant scale; links below open the static portraits. Case numbers identify the saved parameter sets in the linked reports.

## Constructed flows

### Icosahedral flow

![Icosahedral Gaussian flow, cases 177 and 183](docs/assets/portraits/icosahedral.webp)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

[Static portrait](docs/assets/portraits/icosahedral.png)

**Numerical chaos; certified bounds and local equilibria.** A smooth three-state field produces thirty spiral regions around a stable center. Cases 177 and 183 survive 16,000-unit observations, step refinement and independent integration. The 60-edge geometric skeleton has a known integrable counterpart; the sampled return strips are not yet a covering proof.

[Equations and evidence](docs/icosahedral-flow.md) · [Equilibria](docs/icosahedral-structure.md) · [Periodic orbits](docs/icosahedral-cycles.md) · [Return map and prior art](docs/icosahedral-returns.md)

## Numerical algorithms

### RAAR projection algorithm

![Spatial-circle RAAR, cases 1051 and 1483](docs/assets/portraits/raar.webp)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

[Static portrait](docs/assets/portraits/raar.png)

**Certified attracting positive-entropy set.** For case 1051, an explicit strict trapping region contains the certified symbolic dynamics. Its maximal invariant set attracts a neighborhood and has topological entropy at least `log(2)/49` per iteration. A certified fold in this attracting set excludes ambient-neighborhood conjugacy to locally injective maps. Transitivity and equality with the numerical portrait remain unresolved; case 1483 is a numerical companion.

[Trapping theorem and article](docs/raar-trapping.md) · [Fold membership and article](docs/raar-critical-membership.md) · [Symbolic certificate](docs/raar-symbolic-dynamics.md) · [Model and trajectories](docs/raar-circles.md) · [Prior art](docs/raar-priority.md)

### Acceleration memory in a scalar iteration

![Three-history acceleration, cases 54 and 449](docs/assets/portraits/anderson.webp)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

[Static portrait](docs/assets/portraits/anderson.png)

**Numerical chaos with a convergent scalar baseline.** The unaccelerated maps have interval-verified unique fixed points and globally convergent iteration. Adding regularized three-history acceleration produces irregular memory states that persist under 80- and 110-digit checks.

[Model, controls and limitations](docs/anderson-memory.md) · [Reproduction](research/unexpected_search/README.md)

### Adaptive integration of a decaying oscillator

![RK45 controller states, cases 335 and 448](docs/assets/portraits/rk45.webp)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

[Static portrait](docs/assets/portraits/rk45.png)

**Numerical chaos in the solver state.** The physical linear oscillator decays, while RK45's position, velocity and proposed step can remain irregular near the error tolerance. Five-million-step, high-precision and unmodified-library controls agree. RK23 has related retained regimes.

[RK45 results](docs/rk45-oscillator.md) · [RK23 results](docs/rk23-oscillator.md) · [RK45 manuscript](papers/rk45-controller/manuscript.md) · [Prior art](docs/rk45-priority.md)

### Optimizer memory on a scalar quadratic

![Published smooth Yogi variant, cases 2850 and 39](docs/assets/portraits/yogi-tanh.webp)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

[Static portrait](docs/assets/portraits/yogi-tanh.png)

**Numerical evidence for two expanding directions.** The published tanh variant of Yogi retains two positive growth rates in selected regimes under five-million-step, high-precision and TensorFlow checks. Separate results cover standard Yogi, Adamax and AdaBelief. These are regimes of existing algorithms.

[Smooth Yogi](docs/yogi-tanh.md) · [Manuscript](papers/yogi-memory/manuscript.md) · [Standard Yogi](docs/yogi-parameter-plane.md) · [Adamax](docs/adamax-geometry.md) · [AdaBelief](docs/adabelief-quadratic.md)

## Molecular-motor dynamics

### Four-row irregular motion

![Two regimes of the four-row motor model](docs/assets/portraits/motor-four.webp)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

[Static portrait](docs/assets/portraits/motor-four.png)

**Numerical chaos in an eleven-dimensional physical model.** Changing the dimensionless transition-rate asymmetry `b` changes a narrow braid into a broader regime. Long full-spectrum calculations retain a positive leading rate. The portraits show orthonormal displacement contrasts, not three-dimensional embeddings.

[Four-row comparison](docs/motor-four-b-study.md) · [Five-row geometry](docs/motor-five-geometry.md) · [Source equations](docs/motor-rows-source-audit.md)

### Eight-row periodic attraction

![Periodic and irregular eight-row motor trajectories](docs/assets/portraits/motor-eight.webp)

<!-- drawing-date -->
*First recorded portrait: **2026-10-30** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

[Static portrait](docs/assets/portraits/motor-eight.png)

**Certified periodic attraction in 23 dimensions.** At the stated exact parameters, a five-cluster periodic orbit has all 22 transverse Floquet multipliers inside the unit circle. The irregular comparison has six positive finite-time exponents over 20,000 units. The two panels show the same physical projection and scale; only the periodic attractor is certified.

[Periodic-attractor certificate](docs/motor-five-cycle-certificate.md) · [Replication instability](docs/motor-five-replication-limit.md) · [Irregular regimes](docs/motor-activity-continuation.md)

## Research index

The [portrait gallery](docs/attractor-gallery.md) provides static views of the retained visual candidates. This index also includes periodic captures, failed searches, source discrepancies and structural exclusions.

“Certified” refers to the particular interval or analytical statement in a report. Positive finite-time growth and a persistent point cloud are numerical evidence. Neither category establishes historical priority without a separate comparison.

### Constructed flows and maps

| System | Main result | Further evidence |
| --- | --- | --- |
| [Icosahedral Gaussian flow](docs/icosahedral-flow.md) | Numerical chaos; certified absorbing balls | [Equilibria](docs/icosahedral-structure.md), [numerical cycles](docs/icosahedral-cycles.md), [return proposal and known spherical skeleton](docs/icosahedral-returns.md) |
| [Tetrahedral Gaussian flow](docs/tetrahedral-flow.md) | Numerical chaos around a stable center | [Equilibrium certificates](docs/tetrahedral-equilibria.md), [1995 predecessor](docs/tetrahedral-homoclinic-prior.md), [earlier literature audit](docs/polyhedral-priority-followup.md) |
| [Localized Gaussian vortices](docs/gaussian-vortex.md) | Three-state numerical candidate | [Known-shape comparison](docs/attractor-novelty-audit.md) |
| [Dissipative twist map](docs/twist-geometry.md) | Two numerical geometries | [Mechanism and prior art](docs/twist-priority-and-mechanism.md) |

### Numerical algorithms

| System | Main result | Further evidence |
| --- | --- | --- |
| [Spatial-circle RAAR](docs/raar-circles.md) | Certified attracting positive-entropy set at R1051; numerical companion portraits | [Trapping certificate](docs/raar-trapping.md), [fold membership and equivalence obstruction](docs/raar-critical-membership.md), [Symbolic-dynamics certificate](docs/raar-symbolic-dynamics.md), [cycle structure](docs/raar-cycle-structure.md), [folding and preimages](docs/raar-folding.md), [exact structure](docs/raar-structure.md), [prior art](docs/raar-priority.md) |
| [Scalar acceleration memory](docs/anderson-memory.md) | Convergent plain iteration and irregular accelerated memory | [Search and controls](research/unexpected_search/README.md) |
| [Finite-sample FastICA](docs/fastica-dynamics.md) | Four numerical regimes; five-million-step spectra and independent precision controls; one rejected long transient | [Portrait and stopping-rule audit](docs/fastica-dynamics.md), [complete evidence](results/fastica/README.md), [reproduction](research/fastica/README.md) |
| [RK45](docs/rk45-oscillator.md), [RK23](docs/rk23-oscillator.md) | Irregular adaptive-controller states for a decaying oscillator | [RK45 structure](docs/rk45-controller-structure.md), [RK23 structure](docs/rk23-controller-structure.md), [RK45 prior art](docs/rk45-priority.md), [RK23 prior art](docs/rk23-priority.md) |
| [Smooth Yogi](docs/yogi-tanh.md) | Two positive finite growth rates in selected regimes | [Structure](docs/yogi-tanh-structure.md), [manuscript](papers/yogi-memory/manuscript.md) |
| [Yogi parameter plane](docs/yogi-parameter-plane.md) | Retained irregular regimes and regular windows | [Earlier cases](docs/yogi-geometry.md), [smoothing control](docs/yogi-smoothing-control.md), [boundedness theorem](docs/scalar-optimizer-boundedness.md), [prior art](docs/yogi-priority.md) |
| [Adamax](docs/adamax-geometry.md) | Irregular memory and competing cycles | [Structure](docs/adamax-structure.md), [prior art](docs/adamax-priority.md) |
| [AdaBelief](docs/adabelief-quadratic.md) | Numerical scalar-quadratic candidate | [AdaDelta cycle audit](docs/adadelta-cycle-audit.md) |

### Molecular motors

| Result | Evidence |
| --- | --- |
| Exact model and physical normalization | [Motor rows](docs/motor-rows-source-audit.md), [motor pair](docs/motor-pair-source-audit.md) |
| Certified periodic saddle and transverse crossing | [Cycle and least-period audit](docs/motor-cycle-certificate.md), [crossing](docs/motor-crossing-certificate.md), [implementation review](docs/motor-proof-review.md) |
| Certified paired periodic attractor and later saddle | [Paired attractor](docs/motor-paired-certificate.md), [crossing](docs/motor-paired-crossing-certificate.md), [saddle](docs/motor-saddle-certificate.md), [postcritical observations](docs/motor-postcritical.md) |
| Certified five-cluster periodic attraction | [Certificate](docs/motor-five-cycle-certificate.md), [replication limit](docs/motor-five-replication-limit.md), [branch observations](docs/motor-cluster-cycle-branch.md) |
| Numerical irregular geometry | [Four rows](docs/motor-four-geometry.md), [rate-parameter slice](docs/motor-four-b-study.md), [five rows](docs/motor-five-geometry.md), [three asymmetric rows](docs/motor-asymmetry-geometry.md) |
| Collective chaos and hyperchaos candidates | [Activity continuation](docs/motor-activity-continuation.md), [row-count survey](docs/motor-rows-gpu-survey.md), [activity path](docs/motor-rows-activity-path.md), [replication index](docs/motor-replication-index.md) |
| Coordinate and long-time controls | [Motor-pair survey](docs/motor-pair-gpu-survey.md), [physical coordinates](docs/motor-pair-coordinate-audit.md), [log contrasts](docs/motor-log-contrast-audit.md), [long departure](docs/motor-log-seventh-departure.md), [original seed](docs/motor-original-seed-log.md) |
| Additional structural work | [Cluster transverse modes](docs/motor-cluster-transverse.md), [network combinatorics](docs/motor-network-combinatorics.md), [excluded tuple](docs/motor-tuple15-excluded.md) |

### Other physical and biological models

| System | Result and qualifications |
| --- | --- |
| [Synchronematic oscillators](docs/synchronematic-research.md) | Numerical collective motion; [bifurcations](docs/synchronematic-bifurcations.md), [source and priority audit](docs/synchronematic-priority-audit.md) |
| [Elastomer slider](docs/lce-slider-research.md) | Relay-model candidate; [switching-width controls](docs/lce-slider-switching.md), [relay limit](docs/lce-slider-relay-limit.md) |
| [Two-strain epidemic model](docs/two-strain-research.md) | Numerical outbreak irregularity with extremely small prevalence |
| [Coupled optical operators](docs/coupled-operators-research.md) | Irregular histories and periodic captures; [lifetime extension](docs/optical-lifetime-extension.md), [cycle branch](docs/optical-cycle-continuation.md), [arclength continuation](docs/optical-cycle-arclength.md), [saddle departure](docs/optical-saddle-departure.md) |
| Optical model interpretation | [Heat-balance source audit](docs/light-operators-audit.md), [restricted exclusion theorem](docs/coupled-operators-order-audit.md), [known feedback equivalence](docs/coupled-operators-priority-audit.md) |

### Negative, incomplete and excluded searches

These studies are retained to document tested mechanisms and prevent repeated false positives. An incomplete integration is not a negative result.

| System | Outcome or unresolved issue |
| --- | --- |
| [NIFTE thermofluidic engine](docs/nifte-research.md) | Periodic controls and unresolved chaos claims; [physical audit](docs/nifte-physical-priority-audit.md) |
| [Rayleigh spring pair](docs/rayleigh-pair-research.md) | Long transients and weak instability; [structural derivation](docs/rayleigh-pair-structure.md) |
| [HHS mechanochemical cell](docs/hhs-research.md) | Periodic and neutral controls; [source audit](docs/hhs-source-audit.md) |
| [Urease reactor](docs/urease-research.md) | No promoted chaotic attractor |
| [Long-cavity whistle](docs/whistle-screen.md) | No retained chaos among 514 declared trajectories |
| [Ice-cloud microphysics](docs/ice-cloud-screen.md) | Weak signals resolve to numerical cycles |
| [Core–shell gel](docs/core-shell-gel-research.md) | Stiff cases; one endpoint resolves to a [stable periodic orbit](docs/core-shell-gel-stiff-followup.md) |
| [Rigid microfluidic network](docs/microfluidic-research.md) | Published periodic control not recovered; [supplement audit](docs/microfluidic-source-audit.md) |
| [Contact rotors](docs/dcr-research.md) | Continuation ambiguity at strict event boundaries; [source audit](docs/dcr-source-audit.md) |
| [Pulsating heat pipe](docs/heat-pipe-source-audit.md) | Complete source closure unavailable |
| [Solar-array flutter](docs/thermal-flutter-source-audit.md) | Reduced model not uniquely specified by retrieved material |
| [Mineral protocells](docs/mineral-protocell-screening.md) | Planarity or monotone growth excludes the proposed route |
| Further source and structural screens | [Candidate queue](docs/candidate-queue.md), [material memory](docs/material-memory-screening.md), [physical meminductor](docs/meminductor-followup-audit.md), [thermal sheet](docs/thermomechanical-sheet-audit.md), [epithelial chain](docs/epithelial-chain-audit.md) |

### Archive

[All report files](docs) · [Result directories](results/README.md) · [Source provenance](PROVENANCE.md) · [Dependency notices](THIRD_PARTY_NOTICES.md)

Portraits use a common rendering convention and declared coordinates. Dates refer to our drawings, not discovery priority. [Figure provenance](tools/README-portraits.md) · [License](LICENSE)

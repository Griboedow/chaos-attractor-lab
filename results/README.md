# Research evidence

Find the data, checks and proofs behind a portrait. The [main gallery](../README.md) shows the figures and equations. For commands and software, see [reproduction](../docs/reproduction.md). Large trajectories can be recovered separately. [Data storage](../docs/data-storage.md) · [Compact records](compact-data/README.md)

[Recent searches](#recent-shape-searches) · [Posterior maps](#posterior-maps) · [Algorithms](#numerical-algorithms) · [Spatial geometry](#spatial-geometry) · [Molecular motors](#molecular-motors) · [Physical and biological models](#physical-and-biological-models) · [Other data folders](#other-data-folders)

## Recent shape searches

| Date and records | What to look for |
| --- | --- |
| **2026-11-23** | [128,878,020 cycles on one matrix](sinkhorn-choreography/) · [A saddle cycle passes through a sink](cubic-elliptic-sites/resonances.md) · [A straight ray and the surface-wide average](geodesic-typical-growth/) |
| **2026-11-23** | [Chaos on a shared basin boundary](cubic-saddle-passages/) · [A folding loop](cubic-material-lines/) · [Three weak-attraction windows](cubic-elliptic-sites/) · [Mirror spirals in matrix scaling](sinkhorn-route/) · [Motor rows trade roles](motor-symmetry-passages/) |
| **2026-11-22** · [Mathematical art](art_20261122/README.md) | [A seven-lobed exposure shared by two motions](art_20261122/rosette/README.md) · [Five Newton poles](art_20261122/newton/README.md) · [A geodesic mosaic](art_20261122/weave/README.md) · [Seventeen cusps inside the two rings](art_20261122/rings/README.md) |
| **2026-11-22** · [Directional basins and fixed-length instability](attractor_followup_2/README.md) | [Fifteen Newton outcomes on a line](attractor_followup_2/newton/README.md) · [Unbounded expansion per lap](attractor_followup_2/geodesics/README.md) · [Ten phases at an optical saddle](attractor_followup_2/optics/README.md) · [When the plotted orbit enters the rings](attractor_followup_2/annuli/README.md) |
| **2026-11-22** · [Four follow-up proofs](low_hanging_20261122/README.md) | [Binary dynamics inside the rings](low_hanging_20261122/annuli/README.md) · [Exact periods between two circles](low_hanging_20261122/algorithms/README.md) · [Straight and curved geodesics](low_hanging_20261122/geodesics/README.md) · [Three cavity cycles](low_hanging_20261122/optics/README.md) |
| **2026-11-22** · [wow_search_20261122](wow_search_20261122/) | [Four illustrated studies](../docs/chaotic-forms-20261122.md). A single point alternates between [two nested rings](wow_search_20261122/symmetry/README.md) on every step. A [closed straight orbit](wow_search_20261122/lattice/README.md) has a certified transverse multiplier above 4.45 million. Fifteen [Newton basins](wow_search_20261122/boundary/README.md) share one boundary. An [optical-cavity map](wow_search_20261122/orbit/README.md) produces a scalloped fan in a numerical orbit. |
| **2026-11-21** · [shape_search_20261121](shape_search_20261121/) | [Six families of spatial maps and flows](../docs/shape-search-20261121.md). Proofs include [three periodic orbits linked only as a trio](shape_search_20261121/frontiers/cubic/README.md) and a [ten-component attractor containing a binary shift](shape_search_20261121/frontiers/maps/README.md). A proposed [eighth attractor in one flow](shape_search_20261121/frontiers/flows/README.md) remains a numerical candidate. |

## Posterior maps

### Product alphabets

Case 988 has a certified positive-entropy core that is unstable. Whether the observed wings belong to a larger chaotic attractor remains open.

| Records | Finding |
| --- | --- |
| [product-alphabet](product-alphabet/) | Four-wing and multilevel candidates selected from 1,753 inputs. The records include periodic windows and competing fixed points. [Study](../docs/product-alphabet/manuscript.md) |
| [product-symbolic](product-symbolic/) · [product-orientation](product-orientation/) | Coverings give [entropy above 0.0254 per iterate](../docs/product-alphabet/manuscript.md#6-a-certified-positive-entropy-subsystem). A [least-period-92 reversed return](../docs/product-alphabet/manuscript.md#61-a-reversed-return-and-an-equivalence-obstruction) rules out specified ambient conjugacies, including standard three-dimensional Hénon maps and smooth Euclidean flow return maps. |
| [product-mechanism](product-mechanism/) · [product-larger-attractor](product-larger-attractor/) | A certified orbit [escapes the core's covering atlas](product-mechanism/README.md), proving instability. An unstable point has a [validated 200-step continuation near the wings](product-larger-attractor/README.md). Two 50-million-step recurrence runs leave membership in the observed attracting component open. |
| [product-horizon](product-horizon/) · [product-attractor-status](product-attractor-status/) | [Two 100-million-step runs](product-horizon/README.md) include three validated [escapes from the portrait envelope](../docs/product-alphabet/manuscript.md#41-a-hundred-million-step-horizon-and-rare-excursions). All 256 new starts retained positive leading growth over four million updates each. A [singular preimage within 0.000060 of a saved state](product-attractor-status/README.md) excludes the specified trapping ball. Smaller traps remain possible. |
| [product-highorder](product-highorder/) | Nonuniform Cases 6 and 7 retain folded shells over 20 million steps. The regular-grid higher-order cases remain exploratory. [Inputs and portraits](product-highorder/README.md) |

### Circles and spheres

| Records | Finding |
| --- | --- |
| [soft-constellations](soft-constellations/) · [soft-trapping](soft-trapping/) · [soft-symbolic-dynamics](soft-symbolic-dynamics/) | The [Gaussian-posterior circle screen](../docs/soft-constellations.md) led to Case 156. A [trap and entry certificate](../docs/soft-constellation-trapping.md) places the start in an attracting set with [entropy at least log(2)/6](../docs/soft-symbolic-dynamics/manuscript.md). The [comparison with earlier methods](../docs/soft-constellation-priority.md) covers mean shift, MMSE, proximal maps and softmax. |
| [soft-sphere](soft-sphere/) · [canonical-sphere](canonical-sphere/) | The [Gaussian sphere sweep](../docs/soft-sphere/manuscript.md) tested 429 cases. The [canonical posterior search](../docs/canonical-sphere/manuscript.md) tested 656 and retained a folded ribbon and a spatial hyperchaos candidate. |
| [canonical-collision](canonical-collision/) | A trapping box and unstable graph certify [two merged histories and a homoclinic return](../docs/canonical-sphere/manuscript.md#5-certified-attraction-merged-histories-and-a-homoclinic-return) in the maximal attracting set. Its identification with the sampled cloud remains open. |

## Numerical algorithms

| Records | Finding |
| --- | --- |
| [fastica](fastica/) · [fastica-tripling](fastica-tripling/) | A [576-case finite-sample search](../docs/fastica-dynamics.md) includes runs with the FastICA library and higher-precision checks. A separate [constructed FastICA input](../docs/fastica-tripling.md) has a certified trap and positive entropy. |
| [fastica-local-injectivity](fastica-local-injectivity/) · [fastica-preimage-rejection](fastica-preimage-rejection/) | Opposite orientation signs [force failure of local injectivity near the connected attracting set](../docs/fastica-tripling-structure.md). A proposed collision has an [extra preimage outside the attractor](../docs/fastica-preimage-rejection.md), leaving intrinsic equivalence open. |
| [sinkhorn-screen](sinkhorn-screen/) · [sinkhorn-unequal-rates](sinkhorn-unequal-rates/) | In the [576-case scaling screen](../docs/sinkhorn-screen.md), all 40 flagged inputs continued to 5.53 million steps. The [rectangular and unequal-rate study](../docs/sinkhorn-unequal-rates.md) retained two candidates through 106 million steps. [First campaign data](sinkhorn-screen/README.md) · [Unequal-rate data](sinkhorn-unequal-rates/README.md) |
| [adabelief-quadratic](adabelief-quadratic/) | A deterministic optimizer on a quadratic gave positive-growth candidates, a rejected long transient and periodic controls. The original implementation and independent high-precision updates were tested. [Study](../docs/adabelief-quadratic.md) |

## Spatial geometry

| Records | Finding |
| --- | --- |
| [twist-geometry](twist-geometry/) · [twist-followup](twist-followup/) | [Three-axis twist maps](../docs/twist-geometry.md) produced two portraits with eight million iterates each. A [follow-up](../docs/twist-priority-and-mechanism.md) gives an affine embedding of the earlier pumped kicked top and tests its unstable directions. Asymptotic attraction remains open for the numerical candidates. |
| [gaussian-vortex](gaussian-vortex/) · [attractor-novelty](attractor-novelty/) | The [localized-vortex search](../docs/gaussian-vortex.md) includes fresh starts, parameter controls and an adaptive spectrum. The [shape and priority audit](../docs/attractor-novelty-audit.md) compares six candidates with 108 reference clouds. Similarity scores alone establish neither priority nor topological equivalence. |
| [geometry-search-20261031](geometry-search-20261031/) | A [four-link chain candidate](../docs/geometry-search-20261031.md) retained finite-time expansion under independent continuation, but its links frequently crossed. |

## Molecular motors

### Activity and geometry

| Records | Finding |
| --- | --- |
| [motor-rows](motor-rows) · [activity-path](motor-rows/activity-path/) | The [three-row survey](../docs/motor-rows-gpu-survey.md) tested 2,880 settings. The [eight-row activity study](../docs/motor-rows-activity-path.md) retained a numerical hyperchaos candidate. |
| [motor-activity-continuation](motor-activity-continuation/) | [Longer CPU and GPU runs](../docs/motor-activity-continuation.md) include a saddle diagnosis and physical perturbations. The [five-cluster branch](../docs/motor-cluster-cycle-branch.md) loses stability through a complex pair of quotient multipliers near ρ = 2.29823502. Its connection to the irregular regime remains open. [Branch data](motor-activity-continuation/cluster-cycle-branch/) |
| [motor-five-geometry](motor-five-geometry/) · [motor-asymmetry](motor-asymmetry/) | The [five-row ribbon](../docs/motor-five-geometry.md) has cyclic symmetry and numerical phase locking. A separate [asymmetric ribbon](../docs/motor-asymmetry-geometry.md) has physical projections, six spectra and finer runs. Topology remains open. |
| [motor-postcritical](motor-postcritical/) · [motor-pair-coordinates](motor-pair-coordinates/) | [Saddle switching and paired-branch continuation](../docs/motor-postcritical.md) suggested a network. [Direct storage of pair contrasts](../docs/motor-pair-coordinate-audit.md) recovered departures hidden by subtracting nearly equal rows. |
| [motor-log-contrasts](motor-log-contrasts/) · [motor-original-seed-log](motor-original-seed-log/) | [Six departures through 150,000 units](../docs/motor-log-contrast-audit.md) and a [seventh through 280,000](../docs/motor-log-seventh-departure.md) were tracked below ordinary floating-point amplitude ranges. [Seventh-departure manifest](motor-log-contrasts/seventh-review/manifest.json). The [original-seed histories](../docs/motor-original-seed-log.md) disagree at long times, leaving that seed's eventual regime open. |

### Certified cycles, saddles and crossings

| Records | Finding |
| --- | --- |
| [motor-proof](motor-proof/) · [motor-period-proof](motor-period-proof/) · [motor-crossing-proof](motor-crossing-proof/) | The [cycle certificate](../docs/motor-cycle-certificate.md) proves existence and stability properties. Validated orbit tubes fix the least period. A [unique local crossing](../docs/motor-crossing-certificate.md) has exactly six critical directions. |
| [motor-paired-proof](motor-paired-proof/) · [motor-paired-crossing-proof](motor-paired-crossing-proof/) | The [four-pair periodic orbit](../docs/motor-paired-certificate.md) attracts in the full eight-row system. A [local fourfold crossing](../docs/motor-paired-crossing-certificate.md) is certified. The intervening branch and nonlinear bifurcation remain open. |
| [motor-five-cycle-proof](motor-five-cycle-proof/) · [motor-five-replication-proof](motor-five-replication-proof/) | For the [eight-row five-cluster attractor](../docs/motor-five-cycle-certificate.md), all nonflow multipliers have magnitude below 0.13. At the same normalized parameters, [exact copies with 8k rows](../docs/motor-five-replication-limit.md) have at least 2(k−1) unstable directions for k ≥ 2. |
| [motor-saddle-proof](motor-saddle-proof/) · [motor-network-combinatorics](motor-network-combinatorics/) | The [six-cluster saddle](../docs/motor-saddle-certificate.md) has one unstable, one flow and 21 stable directions. The [eight-node graph](../docs/motor-network-combinatorics.md) and its [5,040 labeled saddle images](../docs/motor-saddle-isotropy.md) are exact. Heteroclinic edges remain unproved. [Weighted-volume identity](../docs/motor-volume-contraction.md) |

### Continuum model and source comparisons

- [motor-pde-lifting](motor-pde-lifting/) contains the [density-PDE lifting proof](../docs/motor-pde-lifting.md). Under its assumptions the [compact global attractor](../docs/motor-global-attractor.md) consists of complete bounded first-harmonic trajectories. The [geometry analysis](../docs/motor-global-geometry.md) gives connectedness and dimension lower bounds.
- [motor-priority-followup](motor-priority-followup/) records [source comparisons and cluster-replication corollaries](../docs/motor-followup-priority.md), including 22 polynomial identities.

## Physical and biological models

### Filaments, actuators and optical feedback

| Records | Finding |
| --- | --- |
| [hydrogel-planar-controls](hydrogel-planar-controls/) · [hydrogel-spatial-controls](hydrogel-spatial-controls/) | Six of eight [planar filament scouts](../docs/hydrogel-planar-controls.md) self-intersected. The other two approached numerical cycles. [Spatial controls](../docs/hydrogel-spatial-controls.md) exposed a coarse-mesh equilibrium absent from the continuum model. |
| [hydrogel-slender-screen](hydrogel-slender-screen/) · [hydrogel-variational-checks](hydrogel-variational-checks/) | The [thinner-rod runs](../docs/hydrogel-spatial-controls.md#circular-section-symmetry-and-thinner-filaments) developed subdiameter clearances at stronger fields. In [refinement and continuation](../docs/hydrogel-variational-checks.md), positive finite growth accompanied slow drift. No chaotic attractor was established. |
| [lce-slider](lce-slider/) · [core-shell-gel](core-shell-gel/) | The [liquid-crystal slider](../docs/lce-slider-research.md) gave positive cases outside the source paper's tested stiffness/damping range. [Smoothing the switch](../docs/lce-slider-switching.md) changes the model. An independent [stiff gel continuation](../docs/core-shell-gel-stiff-followup.md) closed on a stable cycle. Other stiff cases remain unclassified. |
| [photofeedback](photofeedback/) · [light-operators](light-operators/) | The [phenomenological feedback model](../docs/photofeedback-research.md) had twelve unsuccessful follow-ups. In a separate optical model the [doubled cycle](../docs/optical-cycle-continuation.md) is unstable. In the [100,000-unit lifetime test](../docs/optical-lifetime-extension.md), both narrow-width cohorts reached cycles and 26 broad-width cases remained finite-time candidates. |

### Fluids, cells and oscillators

| Records | Finding |
| --- | --- |
| [nifte](nifte/) · [rayleigh-pair](rayleigh-pair/) | The [thermofluidic model](../docs/nifte-research.md) recovered stable and unstable source cycles. The [quadratic-drag spring pair](../docs/rayleigh-pair-research.md) gave stable and saddle cycles after 3,081 scouts. No chaos candidate survived. |
| [hhs](hhs/) · [microfluidic](microfluidic/) | None of 186 completed long intervals of the [five-state cell model](../docs/hhs-research.md) retained a chaos candidate. One interval is incomplete. In the [four-state microfluidic model](../docs/microfluidic-research.md) independent solvers failed to recover the source cycle with the printed parameters, so the broader search was deferred. |
| [synchronematic](synchronematic/) · [candidate-research](candidate-research/) | The [phase/orientation study](../docs/synchronematic-research.md) found finite positive growth with the literal article coefficients. The author's code differs in one sign. The quantum-neuron folder includes independent controls and rejected candidates. |
| [dcr](dcr/) · [thermomechanical-sheet](thermomechanical-sheet/) | The [preliminary rotor-network audit](../docs/late-candidate-screening-6.md) identified source/code discrepancies. The [independent DCR controls](../docs/dcr-research.md) distinguish accepted endpoints from callback values. The [sheet audit](../docs/thermomechanical-sheet-audit.md) checks coordinates, Jacobians and equilibria. |

### Screening before simulation

[candidate-audits](candidate-audits/) contains the source equations and the reason each search stopped. The [mineral-membrane analysis](../docs/mineral-protocell-screening.md) gives dimensional and monotonicity exclusions. The [Marangoni models](../docs/late-candidate-screening-5.md) lacked closure or calibration, or already had reported chaos. The [heat-pipe reconstruction](../docs/heat-pipe-source-audit.md) lacked the final model closure. For [thermal flutter](../docs/thermal-flutter-source-audit.md), autonomy and thermal-energy checks left the exact simulation model unspecified.

## Other data folders

The [research index](../README.md#research-index) links the scientific reports. Expand a subject below to open its remaining data folders directly.

<details>
<summary>Projection methods (RAAR)</summary>

[cyclic-raar-search](cyclic-raar-search/) · [raar-circles](raar-circles/) · [raar-constellations](raar-constellations/) · [raar-phase-classes](raar-phase-classes/) · [raar-cycle-structure](raar-cycle-structure/) · [raar-folding](raar-folding/) · [raar-fold-certificate](raar-fold-certificate/) · [raar-critical-membership](raar-critical-membership/) · [raar-trapping](raar-trapping/) · [raar-symbolic-dynamics](raar-symbolic-dynamics/) · [raar-catalog-comparison](raar-catalog-comparison/) · [raar-priority-flow-limit](raar-priority-flow-limit/) · [raar-priority-extension](raar-priority-extension/)

</details>

<details>
<summary>Adaptive integration and optimization</summary>

[rk23-controller](rk23-controller/) · [rk45-controller](rk45-controller/) · [adadelta-cycle-audit](adadelta-cycle-audit/) · [adamax-geometry](adamax-geometry/) · [yogi-geometry](yogi-geometry/) · [yogi-grid-attraction](yogi-grid-attraction/) · [yogi-plane](yogi-plane/) · [yogi-smoothing](yogi-smoothing/) · [yogi-tanh](yogi-tanh/) · [optimizer-priority-followup](optimizer-priority-followup/)

</details>

<details>
<summary>Constructed flows, maps and geometry</summary>

[icosahedral-flow](icosahedral-flow/) · [icosahedral-cycles](icosahedral-cycles/) · [icosahedral-returns](icosahedral-returns/) · [icosahedral-structure](icosahedral-structure/) · [tetrahedral-equilibria](tetrahedral-equilibria/) · [tetrahedral-homoclinic-prior](tetrahedral-homoclinic-prior/) · [polyhedral-priority-followup](polyhedral-priority-followup/) · [coupled-operators](coupled-operators/) · [product-irregular-approach](product-irregular-approach/) · [soft-local-injectivity](soft-local-injectivity/) · [spherical-constellations](spherical-constellations/) · [attractor-geometry](attractor-geometry/) · [unexpected-search](unexpected-search/) · [portrait-presentation](portrait-presentation/)

</details>

<details>
<summary>Molecular motors</summary>

[motor-pair](motor-pair/) · [motor-three-transverse](motor-three-transverse/) · [motor-four-geometry](motor-four-geometry/) · [motor-four-b-study](motor-four-b-study/) · [motor-four-extension](motor-four-extension/) · [motor-ribbon-neighborhood](motor-ribbon-neighborhood/) · [motor-cycle-attraction-visual](motor-cycle-attraction-visual/) · [motor-known-comparison](motor-known-comparison/) · [motor-replication-article](motor-replication-article/)

</details>

<details>
<summary>Physical and biological searches</summary>

[biophysical-screens-20261031](biophysical-screens-20261031/) · [epithelial-chain](epithelial-chain/) · [ice-cloud-screen](ice-cloud-screen/) · [macrospin-exclusion](macrospin-exclusion/) · [meminductor-energy](meminductor-energy/) · [pneumatic-source-audit](pneumatic-source-audit/) · [rheostat](rheostat/) · [spin-crossover](spin-crossover/) · [supramolecular-screen](supramolecular-screen/) · [two-strain](two-strain/) · [urease](urease/) · [whistle-screen](whistle-screen/) · [wing-dsa-screen](wing-dsa-screen/)

</details>

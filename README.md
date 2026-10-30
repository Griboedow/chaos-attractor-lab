# Chaos Attractor Lab

Reproducible investigations of unusual physical and biological oscillator models, using Python, OpenCL parameter surveys and independent CPU checks. **No previously unpublished attractor has been established.** The reports distinguish numerical evidence, known literature results and open hypotheses.

**[Read the findings](docs/research-findings.md)** for the results, negative searches and remaining hypotheses, with figures and primary references.

![Rotating physical projections: certified periodic motor attractor and hyperchaos candidate](docs/assets/motor-regimes-rotation.gif)

Left: the five-cluster periodic attractor. Right: the hyperchaos candidate from the latest CPU continuation. Both use the same physical coordinates and axis limits; the camera rotates around saved trajectories. [Full-size animation](docs/assets/motor-regimes-rotation.gif) · [Still image](docs/assets/motor-regimes-rotation.png) · [Data and reproduction](docs/attractor-geometry.md#motor-regimes-in-physical-coordinates).

![Six nearby full-state trajectories approaching the same periodic motor orbit](docs/assets/motor-cycle-attraction.png)

The dashed black curve is the attracting cycle; the six colored paths start nearby and approach it. The graph measures distance to the orbit in the full state, allowing different phases. [Watch the actual approach](docs/assets/motor-cycle-attraction.gif) · [Full-size picture](docs/assets/motor-cycle-attraction.png) · [Starts, controls and limits](docs/attractor-geometry.md#approach-to-the-periodic-orbit).

## Why this exists

Years ago, a coursework project asked me to scan two parameters of Chua's circuit. I spent months calculating: I understood the subject poorly, and my computer was weak. I returned to the problem with a GPU and better tools to see whether there was something interesting beyond the original assignment. [The Chua atlas](https://github.com/Griboedow/chua-gpu-atlas) keeps that first study; this project collects the subsequent searches in other systems.

## What the calculations found

**A second type of attracting periodic orbit is now certified.** At exact `rho=129/55`, the [five-cluster certificate](docs/motor-five-cycle-certificate.md) proves a nonconstant orbit with multiplicities `2+2+1+2+1`, one neutral flow multiplier and twenty-two multipliers of modulus below 0.13. A changed-precision interval recomputation confirms the result. Its 840 labeled copies are symmetry equivalents of one orbit class; historical novelty and the least return period are not established. A [separate numerical continuation](docs/motor-cluster-cycle-branch.md) finds a complex quotient instability near `rho=2.29823502`, whose relation to the chaotic candidate remains open.

**A local computer-assisted result is now available.** Interval Taylor integration and a Krawczyk inclusion certify a nonconstant four-plus-four periodic motor orbit, its least period and exactly six expanding directions in the full 23-dimensional system. A [stronger crossing certificate](docs/motor-crossing-certificate.md) isolates a unique nearby parameter where exactly six transverse multipliers cross −1 with nonzero speed. The certificates and analytic controls have been recomputed. This proves local properties of periodic motion; it does not prove the numerical hyperchaotic candidate is an attractor. [Base theorem, proof data and reproduction](docs/motor-cycle-certificate.md).

The strongest new numerical lead is **hyperchaos in the published eight-row molecular-motor model**, while varying only its prescribed chemical activity. The independent 23-dimensional CPU calculation now retains six positive exponents over 20,000 observed time units. The complete original GPU cohort retains the same 31 of 32 positives through a further 20,000 units in two numerical settings. A separate audit finds synchronized states and a periodic saddle hidden among the cohort, so positive tangent growth alone does not classify an attractor. The representative lies about 64% above the source's onset activity, beyond its illustrated 10% example. These are finite-time numerical results, with scientific priority unresolved. [Long continuation, physical perturbations and limitations](docs/motor-activity-continuation.md), [original activity-path study](docs/motor-rows-activity-path.md).

A second [interval certificate proves an attracting four-pair periodic orbit](docs/motor-paired-certificate.md): one neutral flow direction and twenty-two contracting directions in the full eight-row model. A [further local certificate](docs/motor-paired-crossing-certificate.md) proves a unique fourfold crossing through minus one near `rho=3.2748747166`, with nonzero speed and certified stability counts on both sides. The [numerical continuation](docs/motor-postcritical.md) connects the observed regimes, but a rigorous connecting branch and nonlinear bifurcation theorem remain open. Beyond that, apparent regular states are periodic saddles, with evidence for switching between symmetry-related copies. Rounding tiny row differences to zero can artificially stop these switches and hide an expanding tangent direction; the report includes a controlled demonstration.

An [exact change of all 23 coordinates](docs/motor-pair-coordinate-audit.md) retains those tiny differences and resolves four successive departures over 40,000 time units. Two integration settings agree on the switches, including recovery from a difference near `10^-75`; a later difference near `10^-140` remains nonzero and regrows. The [network proof strategy](docs/motor-network-proof-strategy.md) separates this finite evidence from the missing connection and attraction theorems, and records established switching precedents.

The [logarithmic continuation](docs/motor-log-seventh-departure.md) now reproduces seven departures through 280,000 time units in two integration settings. A contrast near `10^-486` recovers and separates its pair again; another reaches `10^-911` while remaining represented. Increasing residence times approach the ratio predicted from the certified local saddle rates. The earlier [150,000-unit study](docs/motor-log-contrast-audit.md) remains frozen, and the global connection and attraction assumptions remain open.

The [six-cluster saddle itself is now certified](docs/motor-saddle-certificate.md), including its least period and the full count of one expanding, one flow and twenty-one contracting directions. An [exact finite-group calculation](docs/motor-network-combinatorics.md) identifies eight distinct symmetry images in the group generated by the observed departures, with a sixteen-edge candidate graph. The periodic states are certified; the connecting trajectories and attraction of the network remain open.

There is also a concrete mechanism clue. The source's periodic four-plus-four cluster motion loses stability to row differences at about 24% above onset. Six transverse directions become unstable, while the periodic motion stays stable inside its cluster subspace. The matching six expanding directions of the later hyperchaotic candidate motivate a cluster-breakup hypothesis; the global connection is still unproved. [Periodic-cluster continuation and independent checks](docs/motor-cluster-transverse.md).

The **three-row version** provides a striking physical-coordinate projection. Three starts retain positive expansion through another 10,000 units, with every one of 80 windows per start positive. A separate CPU calculation gives one positive, one nearly neutral and six negative exponents. The rows do not collapse into a fixed pair. [Full evidence, failed candidates and reproduction](docs/motor-rows-gpu-survey.md).

![Rotating physical projection of an eight-dimensional motor-row candidate](results/motor-rows/cpu-candidate-N3-783.gif)

The curved triangular envelope is a projection of eight state variables onto three physical observables. Permutation symmetry helps explain its appearance; an unusual picture alone does not establish a new attractor class.

The [follow-up priority and symmetry audit](docs/motor-followup-priority.md) derives exact replication and cluster-spectrum identities, including the target periodic saddle's `N−2` expanding directions for every even row count. It also records established cluster-breakup precedents and two unread author theses. These algebraic consequences and source gaps remain separate from an exact-model novelty claim.

The [density-PDE lifting argument](docs/motor-pde-lifting.md) shows why the specified sinusoidal source model closes exactly: omitted density modes contract exponentially and add only stable directions to a periodic orbit. The local cycle theorems therefore extend to that density PDE under its stated assumptions. This does not extend them to a bending flagellum coupled to fluid or to arbitrary transition-rate laws.

The [global dissipativity theorem](docs/motor-global-attractor.md) also proves a compact global attractor and identifies the PDE attractor exactly with the finite first-harmonic lift. Its existence does not determine whether the attracting dynamics are periodic, chaotic or a switching network.

An [exact weighted-volume identity](docs/motor-volume-contraction.md) further proves zero phase-space volume of that global attractor wherever its stated coefficient is negative, including the main hyperchaotic candidate and the certified switching saddle. This sufficient condition does not cover the lower-rho lane-46 and five-cluster-cycle parameters. The identity bounds the sum of expansion rates while allowing several individual positive rates.

The clearest lesson is how easily a long chaotic transient can resemble an attractor. In the optical model, one trajectory kept positive finite-time expansion for 72,875 thermal-time units before reaching a checked periodic cycle. At both narrow optical widths, all 32 initial conditions eventually settled into cycles. At the broad width, 26 of 32 remained irregular through 100,000 units; the other six had already reached cycles. These are observations of a fixed cohort under an assumed logistic illumination law, not an experimental calibration or proof of an asymptotic chaotic attractor. [Read the long-run report](docs/optical-lifetime-extension.md).

![Long optical transients and persistence through 100,000 thermal-time units](docs/assets/optical_lifetime_extension.png)

The periodic destination also has structure: continuation finds a branch that loses and regains stability, together with an unstable doubled branch. [The bifurcation report](docs/optical-cycle-continuation.md) separates those local observations from claims about global basins or chaotic attractors.

A [further continuation](docs/optical-cycle-arclength.md) locates a fold at optical width 0.0404788044. Stable and saddle cycles coexist at width 0.04; both are independently checked. Their connection to the basin boundary of the irregular regime remains open.

A [two-sided departure test](docs/optical-saddle-departure.md) now gives a concrete clue: small perturbations near that saddle lead toward the stable cycle on one side and retain irregular motion on the other, at both tested perturbation sizes. Independent Radau runs confirm the cycle-proximity outcome for one size; the second-size audit is explicitly incomplete. This is finite evidence for a basin-boundary hypothesis, not a proof of a new attractor.

![Physical projection of two optical departure outcomes](docs/assets/optical_saddle_departure_rotation.gif)

## Start with the evidence

| Status | Models and reports | Interpretation |
| --- | --- | --- |
| Certified local periodic dynamics | [Motor-cycle certificate](docs/motor-cycle-certificate.md), [unique crossing](docs/motor-crossing-certificate.md), [implementation review](docs/motor-proof-review.md) | Validated periodic orbit, least period, six expanding directions, and a unique local transverse crossing with nonzero speed. The global route to chaos and scientific priority remain open. |
| Attracting periodic motion and later saddle switching | [Four-pair certificate](docs/motor-paired-certificate.md), [saddle certificate](docs/motor-saddle-certificate.md), [postcritical study](docs/motor-postcritical.md) | A locally attracting cycle and a later periodic saddle are certified in all 23 dimensions. Repeated switching has numerical evidence; infinite-time connections and network attraction are unproved. |
| Attracting five-cluster periodic motion | [Five-cluster certificate](docs/motor-five-cycle-certificate.md), [replication limit](docs/motor-five-replication-limit.md), [numerical branch](docs/motor-cluster-cycle-branch.md) | Full 23D attraction is certified at exact rho129/55. Exact copies with N=8k, k>=2, are transversely unstable at the same normalized parameters. A lower-rho complex quotient instability is numerical; torus birth, connection to chaos and historical novelty remain open. |
| Molecular-motor chaos candidates | [Motor-pair survey](docs/motor-pair-gpu-survey.md), [exact source audit](docs/motor-pair-source-audit.md) | Eight trajectories retain positive estimates through tighter 24,063-unit observations; an independent CPU run supports positive growth. The neutral exponent still needs convergence. The model already appeared in 2014; scientific priority and attracting-chaos claims remain open. |
| Collective motor chaos and hyperchaos candidates | [Long activity continuation](docs/motor-activity-continuation.md), [broader N-row survey](docs/motor-rows-gpu-survey.md) | Eight rows: six positive directions through 20,000 observed CPU units; 31/32 original GPU histories retain positive growth through a further 20,000 units in two settings, with synchronized cases separately audited. Three rows: one positive direction and 10,000-unit persistence. Scientific priority and asymptotic attraction remain open. |
| Numerical chaos candidates requiring model/priority qualifications | [Synchronematic oscillators](docs/synchronematic-research.md), [LCE slider](docs/lce-slider-research.md), [two-strain model](docs/two-strain-research.md) | Read each report's source variants, physical range, solver checks and prior-art limits. A positive finite-time exponent alone is insufficient. |
| Coupled optical operators: persistent candidates and long chaotic transients | [Two light-mediated material operators](docs/coupled-operators-research.md) | The continuous heat balance of SI Eq. (10) keeps temperature continuous at switches; the author's implementation was not recovered. Smooth optical profiles are explicitly assumed extensions. Both narrow-width cohorts reach independently checked periodic cycles; 26 broad-width trajectories remain irregular through 100,000 thermal-time units. These counts are finite-time cohort observations. |
| Negative or unresolved searches | [NIFTE thermofluidic engine](docs/nifte-research.md), [Rayleigh spring pair](docs/rayleigh-pair-research.md), [HHS mechanochemical cell](docs/hhs-research.md), [urease reactor](docs/urease-research.md) | Long transients, weak instability and incomplete integrations are recorded explicitly; no retained chaos is promoted to an attractor. |
| Resolved stiff endpoint | [Core–shell gel follow-up](docs/core-shell-gel-stiff-followup.md) | One saved step-budget endpoint closes to a strongly stable periodic orbit under independent implicit CPU solvers; the other unresolved cases remain unclassified. |
| Structural exclusions and source audits | [Candidate queue](docs/candidate-queue.md), [material-memory screening](docs/material-memory-screening.md), [physical meminductor](docs/meminductor-followup-audit.md) | Mathematical reductions and missing physical closure can rule out or postpone a search. |
| Source-control discrepancy | [Rigid microfluidic network](docs/microfluidic-research.md), [final supplement audit](docs/microfluidic-source-audit.md) | The printed closure approaches a stable equilibrium in our stated controls; the source's periodic example was not recovered. No coefficient was refitted and no mass chaos scan was started. |
| Event law unresolved | [Motor-contact rotors](docs/dcr-research.md), [contact-source audit](docs/dcr-source-audit.md) | Short independent implicit controls are archived, but strict contact boundaries leave a continuation ambiguity. Budget-limited solver callbacks are not accepted trajectory endpoints or attractor evidence. |
| Restricted structural exclusion | [Optical-order audit](docs/coupled-operators-order-audit.md) | For the explicitly smooth optical law and real mechanical poles, the applicable theorem excludes equilibrium-free chaotic limit sets. Its scope does not cover the underdamped candidates or ideal relay. |
| Priority remains unresolved | [Optical relay and inertial-neuron comparison](docs/coupled-operators-priority-audit.md) | An exact tanh normalization connects this model to established feedback classes. The bounded literature comparison does not establish a new mechanism, topological attractor class or exact-model priority. |
| Periodic-branch mechanism | [Optical-cycle continuation](docs/optical-cycle-continuation.md) | A periodic branch loses and regains stability, with a numerically unstable doubled branch. These local bifurcation observations do not establish an attracting chaotic set. |
| Source-only screening | [Pulsating heat pipe](docs/heat-pipe-source-audit.md) | The complete 2026 closure was not retrieved; the accessible 2024 predecessor already discusses prospective chaotic motion. No reconstructed 2026 model was simulated. |
| Source-definition gate | [Solar-array thermal flutter](docs/thermal-flutter-source-audit.md) | An unambiguous autonomous reduced model and full source-control settings were not recovered. A thermal-subsystem energy identity is restricted to its stated physical domain. |
| Structural exclusion | [Mineral protocells](docs/mineral-protocell-screening.md) | The latest model is planar; its growing-membrane predecessor has monotone thickness in the active interior. No new trajectory search was run. |

![Synchronematic evidence](docs/assets/synchronematic_evidence.png)

This figure is a projection of a higher-dimensional state, not a claim of a new three-dimensional attractor. [Geometry and symmetry notes](docs/attractor-geometry.md) explain that distinction.

![NIFTE controls](docs/assets/nifte-controls.png)

The NIFTE controls include stable and unstable periodic orbits. An unstable periodic orbit is not an attracting chaotic set.

![Rotation of physical observables](docs/assets/synchronematic_projection_rotation.gif)

The full source Fig. 37 periodic control, independent CPU spectra, optical-width sensitivity and 54 original capture checks remain in [the operator report](docs/coupled-operators-research.md). The [extension](docs/optical-lifetime-extension.md) adds 16 independently checked captures while preserving original histories and settings. The 26 surviving trajectories are right-censored; these counts are not basin probabilities or an escape-rate fit.

## Run locally

The later motor modules were also checked in the updated independent public clone at `ad83bf2`: N=3/N=8 GPU versus CPU controls passed under isolated Python, with no Chua imports or changes to archived results. [Additional packaging record](provenance/motor-public-clone-validation.json).

The published five-cluster proof tree at `bcd5119` was separately exported from the fetched public branch. All 24 manifest hashes matched; its algebra controls and exact symmetry enumeration passed in the separate proof environment. Local links were checked, while browser rendering was unavailable for these latest reports. [Publication verification](provenance/motor-session-publication-20261030.json).

**Public-clone check:** a fresh clone of this repository at `a514509`, with a newly installed environment, passed four GPU/CPU runtime controls, the optional DCR control and `pip check`. The checks ran outside both repositories with isolated Python, imported no Chua module and left the archived evidence unchanged. [Validation record](provenance/public-clone-validation.json), [file-copy audit](provenance/publication-copy-integrity.json). This verifies the tested packaging paths, not every scientific result or platform.

From the project root, create an environment with `python -m venv .venv`. Activate it using `source .venv/bin/activate` on Unix-like systems or `.venv\Scripts\Activate.ps1` in PowerShell, then run `python -m pip install -r requirements.txt`. Python 3.14.3 was used for the recorded work. A separate fresh environment passed installation, `pip check` and four extracted GPU/CPU controls at source snapshot `8a3564c`; those tested modules and the runtime extraction are unchanged in this refreshed copy. The [validation record](provenance/packaging-validation.json) retains those checks, and [requirements-lock.txt](requirements-lock.txt) records the observed environment. This does not certify every model, every platform or a future public clone. GPU scripts require a working OpenCL GPU driver with FP64 support. They do not silently substitute CPU execution.

The optional archived DCR source-solver controls additionally require `python -m pip install -r requirements-dcr.txt` after installing the main requirements. This adds pinned `scikit-sundae`/SUNDIALS bindings. Version 1.2.0 passed a separate fresh-environment dependency check and an isolated CPU control of duration 0.01. Its maximum state discrepancy from the archived fine control was 1.16e-10; [optional-dependency provenance](provenance/optional-dcr-packaging.json) retains the original 3cf58d0 snapshot and confirms the relevant files are unchanged here. This is a packaging check, not a new scientific result. The observed lock includes this optional package. It was not part of the four-module base-runtime validation and no DCR mass survey is supported. See [the DCR reproduction instructions](docs/dcr-research.md) before starting its bounded controls.

Run existing command-line scripts from the repository root: some default result paths depend on the current working directory. Examples using saved data, without a new parameter survey:

```sh
python research/synchronematic_observables.py
python research/plot_synchronematic.py
python research/plot_nifte.py
python research/plot_rayleigh_pair.py
python research/plot_hhs.py
python research/plot_microfluidic.py
python research/plot_core_shell_stiff.py
python research/extend_optical_lifetimes.py --plot
```

Plotting commands write figures. Keep the archived evidence in a clean checkout and use a disposable copy when reproducing figures or simulations. The [results index](results/README.md) and individual reports identify settings, seeds, traces and exact model versions. Many simulation scripts overwrite named result files; read their arguments before running.

## Provenance and limitations

This project was extracted from a fixed [Chua atlas snapshot](https://github.com/Griboedow/chua-gpu-atlas/tree/b60418fc2fbdd8a34ab873c5581e31f1b9c0727a). See [provenance](PROVENANCE.md), [MIT license](LICENSE) and [dependency notices](THIRD_PARTY_NOTICES.md). The original project's historical commits remain intact. Existing rendering records certify their stated old-repository revisions; independent checks of this public repository will be recorded separately. The verified forward transfer is complete. Each later science block records its own source provenance.

## References

1. Guo, H., Li, K., Yang, J., Li, D., Liu, F., and Zeng, H. (2025). Light-mediated communication in responsive materials ranging from individual self-oscillators to feedback-driven network. *Nature Communications* **16**, 11315. [doi:10.1038/s41467-025-66395-3](https://doi.org/10.1038/s41467-025-66395-3). Physical source; the logistic optical profile used here is an explicit modeling assumption.
2. Leyva, S. G., Zhang, Z., Olvera de la Cruz, M., and Bishop, K. J. M. (2026). Self-oscillating synchronematic colloids. *Nature Communications* **17**, 1841. [doi:10.1038/s41467-026-68552-8](https://doi.org/10.1038/s41467-026-68552-8). Source coefficient variants are audited in the accompanying reports.
3. Wang, Y., Markides, C. N., and Chachuat, B. (2019). Optimization-based investigations of a two-phase thermofluidic oscillator for low-grade heat conversion. *BMC Chemical Engineering* **1**, 20. [doi:10.1186/s42480-019-0019-3](https://doi.org/10.1186/s42480-019-0019-3). Source of the NIFTE equations and periodic controls.
4. Letellier, C., and Aguirre, L. A. (2012). Required criteria for recognizing new types of chaos: Application to the “cord” attractor. *Physical Review E* **85**, 036204. [doi:10.1103/PhysRevE.85.036204](https://doi.org/10.1103/PhysRevE.85.036204). Motivation for distinguishing an unusual projection from a new dynamical class; no such topological classification is established here.
5. Alouges, F., Anello, I., DeSimone, A., Lefebvre-Lepot, A., and Levillain, J. (2025). Some mathematical models for flagellar activation mechanisms. *Mathematical Models and Methods in Applied Sciences* **35**, 2395–2424. [doi:10.1142/S0218202525500423](https://doi.org/10.1142/S0218202525500423). Source of the exact constrained motor-row equations.

Primary references and precise equation variants accompany each system report, including [synchronematic source/priority checks](docs/synchronematic-priority-audit.md), [optical heat-balance and switching audit](docs/light-operators-audit.md), [NIFTE physical and source audit](docs/nifte-physical-priority-audit.md), [Rayleigh structural derivation](docs/rayleigh-pair-structure.md), [HHS source audit](docs/hhs-source-audit.md) and [microfluidic final-source check](docs/microfluidic-source-audit.md). The [candidate queue](docs/candidate-queue.md) provides the wider literature trail. Referencing a paper does not establish priority for our observations; failure to find a publication is not evidence that none exists.

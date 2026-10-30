# Chaos Attractor Lab

Reproducible investigations of unusual physical and biological oscillator models, using Python, OpenCL parameter surveys and independent CPU checks. **No previously unpublished attractor has been established.** The reports distinguish numerical evidence, known literature results and open hypotheses.

**[Read the findings](docs/research-findings.md)** for the results, negative searches and remaining hypotheses, with figures and primary references.

## Why this exists

Years ago, a coursework project asked me to scan two parameters of Chua's circuit. I spent months calculating: I understood the subject poorly, and my computer was weak. I returned to the problem with a GPU and better tools to see whether there was something interesting beyond the original assignment. [The Chua atlas](https://github.com/Griboedow/chua-gpu-atlas) keeps that first study; this project collects the subsequent searches in other systems.

## What the calculations found

**A local computer-assisted result is now available.** Interval Taylor integration and a Krawczyk inclusion certify a nonconstant four-plus-four periodic motor orbit, its least period and exactly six expanding directions in the full 23-dimensional system. A [stronger crossing certificate](docs/motor-crossing-certificate.md) isolates a unique nearby parameter where exactly six transverse multipliers cross −1 with nonzero speed. The certificates and analytic controls have been recomputed. This proves local properties of periodic motion; it does not prove the numerical hyperchaotic candidate is an attractor. [Base theorem, proof data and reproduction](docs/motor-cycle-certificate.md).

The strongest new numerical lead is **hyperchaos in the published eight-row molecular-motor model**, while varying only its prescribed chemical activity. An independent 23-dimensional CPU calculation retains six substantial positive exponents over 4,000 observed time units. Separately, 31 of 32 selected GPU trajectories retain positive growth through another 5,000 units. The representative lies about 64% above the source's onset activity, beyond its illustrated 10% example. These are finite-time numerical results, with scientific priority unresolved. [Activity-path evidence and limitations](docs/motor-rows-activity-path.md).

There is also a concrete mechanism clue. The source's periodic four-plus-four cluster motion loses stability to row differences at about 24% above onset. Six transverse directions become unstable, while the periodic motion stays stable inside its cluster subspace. The matching six expanding directions of the later hyperchaotic candidate motivate a cluster-breakup hypothesis; the global connection is still unproved. [Periodic-cluster continuation and independent checks](docs/motor-cluster-transverse.md).

The **three-row version** provides a striking physical-coordinate projection. Three starts retain positive expansion through another 10,000 units, with every one of 80 windows per start positive. A separate CPU calculation gives one positive, one nearly neutral and six negative exponents. The rows do not collapse into a fixed pair. [Full evidence, failed candidates and reproduction](docs/motor-rows-gpu-survey.md).

![Rotating physical projection of an eight-dimensional motor-row candidate](results/motor-rows/cpu-candidate-N3-783.gif)

The curved triangular envelope is a projection of eight state variables onto three physical observables. Permutation symmetry helps explain its appearance; an unusual picture alone does not establish a new attractor class.

The [follow-up priority and symmetry audit](docs/motor-followup-priority.md) derives exact replication and cluster-spectrum identities, including the target periodic saddle's `N−2` expanding directions for every even row count. It also records established cluster-breakup precedents and two unread author theses. These algebraic consequences and source gaps remain separate from an exact-model novelty claim.

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
| Molecular-motor chaos candidates | [Motor-pair survey](docs/motor-pair-gpu-survey.md), [exact source audit](docs/motor-pair-source-audit.md) | Eight trajectories retain positive estimates through tighter 24,063-unit observations; an independent CPU run supports positive growth. The neutral exponent still needs convergence. The model already appeared in 2014; scientific priority and attracting-chaos claims remain open. |
| Collective motor chaos and hyperchaos candidates | [Published activity path](docs/motor-rows-activity-path.md), [broader N-row survey](docs/motor-rows-gpu-survey.md) | Eight rows: six positive directions in a full 23D CPU spectrum; 31/32 selected trajectories persist for 5,000 additional units. Three rows: one positive direction and 10,000-unit persistence. Scientific priority and asymptotic attraction remain open. |
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

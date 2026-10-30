# Priority and positioning

Audit date: 30 October 2026. Revision 2 expands the comparison to equivalent terminology and separates inherited theory from the contribution proved here. This is a bounded literature audit, not a certificate of historical priority.

## What the paper contributes

**An explicit, interval-certified periodic orbit and the complete stability classification of its balanced replicas in an existing molecular-motor model.** At specified parameters, the eight-row orbit has five distinct groups of sizes `2+2+1+2+1` and is attracting in the full 23-dimensional state space. For every integer `k >= 1`, the same weighted motion at `N=8k` has the exact Floquet index

```text
(unstable, neutral, stable) = (2k-2, 1, 22k).
```

All stable multipliers have modulus below `0.13`. Every unstable multiplier is real and lies in `(-23.084,-23.083)`. The contribution consists of the enclosed orbit, the certified block spectra, their complete count for an unbounded discrete family, and the reproducible proof packet. The qualitative parameter-persistence corollary follows by standard arguments; it is not presented as a new persistence theorem.

**The general mechanism is established background.** The paper is not a discovery of network enlargement, cluster splitting, absent singleton modes, or transverse instability. Its explicit decomposition is retained as a short model-specific derivation needed to connect the equations to the certificate. No repeated numerical study of the general phenomenon is needed.

## Equivalent terms and closest precedents

| Our terminology | Established terminology and primary source | Consequence for positioning |
| --- | --- | --- |
| Balanced replication preserving group motion | Network inflation; quotient lifting. [Aguiar, Ashwin, Dias and Field, 2011](https://doi.org/10.1007/s00332-010-9083-9), Section 3.4; [author manuscript](https://cmup.fc.up.pt/cmup/apdias/cycles.pdf) | Retaining dynamics in a synchrony subspace of a larger network is known. Our weighted mean-field construction is a specific realization. |
| Internal differences; no splitting direction in a singleton | Tangential/transverse cluster stability. [Orosz, Moehlis and Ashwin, 2009](https://doi.org/10.1143/PTP.122.611), Eqs. (2.19)-(2.23), pp. 615-616; [published paper](https://sites.me.ucsb.edu/~moehlis/moehlis_papers/PTPpublished.pdf) | The multiplicity `cluster size - 1`, including zero for singletons, is explicit prior work. We cite this reasoning rather than claim it. |
| Stability of a pair contrast | Evaporation or split Lyapunov exponent. [Zaks and Pikovsky, 2017](https://doi.org/10.1038/s41598-017-04918-9), Methods; [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5498578/) | Opposite perturbations of synchronized units are an established stability diagnostic. |
| Quotient and transverse blocks | Symmetry-based cluster synchronization and isolated desynchronization. [Pecora et al., 2014](https://doi.org/10.1038/ncomms5079); [Sorrentino et al., 2016](https://doi.org/10.1126/sciadv.1501737) | The framework is inherited. The particular motor matrices and strict bounds must still be established. |
| Apparent irreversible numerical synchronization | Finite-precision merging of identical units. [Zaks and Pikovsky, 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5498578/), Methods, “Integration of equations” | The numerical pitfall in the separate switching programme is already described. Recovery of motor contrasts below ordinary floating-point range is an application, not a new-pitfall claim. |

For a periodic internal block, an evaporation exponent is `log(max |mu|)/T`. It identifies growth but discards multiplier sign and does not count every stable and unstable direction. Our complete three-dimensional spectra provide both. A real scalar internal phase equation has multiplier `exp(integral a(t) dt) > 0`; that scalar approximation cannot reproduce our negative internal multiplier. This distinguishes the calculations without suggesting that negative multipliers were previously unknown.

Balanced replication keeps `m_a/N` fixed and changes block multiplicities from `m_a-1` to `k*m_a-1`. Once the block spectra are known, the index formula is dimension counting. The substantive work is certifying the orbit and those spectra in the fixed motor equations; the counting principle itself is not a discovery.

## Exact model equivalence versus related phenomena

| Comparison | Equivalence level | What is inherited / what remains specific |
| --- | --- | --- |
| [Alouges et al., M3AS 2025](https://doi.org/10.1142/S0218202525500423), [full preprint v2](https://arxiv.org/html/2409.03506v2) | Same first-harmonic motor equations after normalization | The model, closure and published two-group numerical examples are inherited. The inspected article does not supply this five-cluster orbit and complete interval-certified replication index. |
| [Oriola et al., EPL 2014](https://doi.org/10.1209/0295-5075/107/18002), Eq. (10) | Exact invertible rescaling of the closed **two-row three-state core** | Earlier periodic/subharmonic motor motion must be credited. The [coordinate map and checks](../../docs/motor-pair-source-audit.md) are already recorded; new notation does not create a new model. This map does not identify the full 14-dimensional five-cluster quotient with that core. |
| [Haugland and Krischer, PRE 2021](https://doi.org/10.1103/PhysRevE.103.L060201) | Related cluster symmetry breaking in Stuart-Landau equations | Related mechanisms do not establish exact equation equivalence. No conjugacy to the full motor vector field has been established. We claim neither a new generic route to chaos nor equivalence on visual resemblance alone. |
| [Kori and Kuramoto, PRE 2001](https://doi.org/10.1103/PhysRevE.63.046214) | Prior slow switching in coupled oscillators | Background for the separate switching programme. Our unproved connecting trajectories cannot be promoted to a heteroclinic-network theorem. |

Terms such as cluster integrity, conditional stability, blow-up, equitable partition and graph fibration can help discover related constructions; they are not additional contributions to claim. Exact equation equivalence requires an explicit map, parameter correspondence and time transformation. A local normal form or a similar projected curve is insufficient. No exhaustive exclusion of all conjugate formulations has been established.

## Defensible novelty statement

The inspected sources do not identify the same explicit motor orbit together with the complete certified index and strict spectral bounds. This supports positioning the paper as a **specific computer-assisted result in the published model**. It does not establish that nobody has computed a similar orbit, or that every component of the proof is new.

> Using established cluster-stability theory, we establish an explicit computer-assisted periodic-orbit and replication-stability result in the first-harmonic motor model: a fully attracting eight-row five-cluster cycle whose balanced replicas have the exact Floquet index `(2k-2,1,22k)` for all `N=8k`.

Avoid “a new general replication mechanism,” “the first cluster desynchronization,” “the first motor chaos,” “a new synchronization-roundoff failure,” “840 new attractors,” or “a proved strange attractor.” The theorem proves a periodic attractor at `N=8` and periodic saddles at larger replicated sizes. Its displayed period is a certified positive return period, not a proved least period.

## Remaining priority limits

- Irene Anello, *Mathematical modeling of spontaneous oscillations in flagellar and bio-inspired systems* (2025), [institutional record](https://iris.sissa.it/handle/20.500.11767/147991). The relevant full chapters remain unread: ordinary institutional PDF retrieval previously returned HTTP 403. Indexed text and metadata are not a substitute.
- Jessie Levillain, *Filament models for swimming at the microscopic scale* (2024), [official HAL record](https://theses.hal.science/tel-04766775v1), [author page](https://jesnial.github.io/). The relevant full chapter remains unread after the earlier access challenge.

These are identified gaps, not an exhaustive list of possible prior sources. Stronger “first” wording is not justified by the absence of an exact phrase match. No author has been contacted. The [search record](literature-search-log.md) states the scope and evidence locations so future work can focus on the remaining comparisons instead of repeating established theory.

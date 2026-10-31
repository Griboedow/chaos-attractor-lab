# Biophysical searches that did not retain a chaos candidate

The completed homogeneous Rho–actin and microbial cross-feeding surveys do not supply a validated strange-attractor candidate. Weak positive finite-time rates mostly accompany ordinary numerical periodic returns. A motor–clutch model remains at a source-interpretation gate. These records preserve the investigated parameter ranges, unresolved cases and numerical failures so that later searches do not mistake them for new leads or repeat them silently.

The [frozen packet](../results/biophysical-screens-20261031/) contains complete scout arrays, declared inputs, executed source snapshots, independent controls, return diagnostics and a hash inventory of the larger local histories. Rechecks overlap earlier designs; their trajectory counts must not be added as independent parameter coverage.

## Rho activation, substrate depletion and actin inhibition

The source is [Chomchai et al., Current Biology (2025)](https://doi.org/10.1016/j.cub.2025.02.003), Supplementary Equation 7 and Table S2, checked in the [author-hosted full paper](https://biology.ed.ac.uk/sites/default/files/2025-04/70-Goryachev-Fcurbiol0325.pdf). We implement its homogeneous three-state kinetics, without diffusion:

```math
\dot R_T=a(R_T)R_D-b(F)R_T,\qquad
\dot R_D=-a(R_T)R_D+b(F)R_T+k_5-k_6R_D,
```

```math
\dot F=k_7+\frac{k_8R_T^2}{1+k_9R_T^2}-k_{10}F,\quad
a=k_0+\frac{\mathrm{GEF}\,k_1R_T^3}{1+k_2R_T^2},\quad b=k_3+k_4\mathrm{GAP}\,F.
```

The printed `k6=.0208333` is retained. The computational time is `tau=.025 t`, where the source time `t` is in seconds. The source homogeneous control approaches an equilibrium near `(2.1513875,3.0000048,12.142266)`; the source paper's spatial waves are a different experiment. Positive concentrations are structurally preserved, but we have not proved global boundedness over the searched parameter family.

| Survey | Declared lanes | Completed | Unresolved step budgets | Primary retained |
| --- | ---: | ---: | ---: | ---: |
| Four-parameter pilot | 2,050 | 1,715 | 335 | 0 |
| Focused unstable-equilibrium screen | 12,860 | 12,756 | 104 | 0 |
| Expanded source-kinetics screen | 8,016 | 7,883 | 133 | 0 |

The pilot varies GEF, GAP, actin speed and exchange speed. The focused study selects 6,430 of 16,384 proposed tuples using a declared equilibrium heuristic, then uses two starts per tuple. Its warm-up and observation are 1,000 and 4,000 tau. The expanded study selects 4,008 of 32,768 proposed tuples while additionally varying basal activation, inactive-pool supply and both saturation constants; its warm-up and observation are 500 and 2,000 tau. Exact logarithmic bounds, seeds, initial states and filters are in each archived protocol. The equilibrium filters are search heuristics, not exclusions of other attractors.

Primary retention requires all four leading-rate windows above `.005` per tau. An exploratory `.001` threshold leaves 32 focused and 97 expanded weak cases. In every expanded weak case, the final tangent is aligned with the autonomous flow to within `2.81e-12` in `1-|cos(angle)|`. That endpoint diagnostic is consistent with phase growth and does not itself determine an asymptotic exponent.

All 32 focused weak cases receive an independent 1,000-tau Radau return screen. **Thirty have lag-one numerical returns; two, IDs 1238 and 4166, provide only four events and remain inconclusive.** Four selected expanded-kinetics weak cases also have lag-one returns, with maximum scaled residual below `5.76e-9`. These are numerical recurrence checks, not interval orbit certificates.

Two pilot weak cases receive 8,000-tau full-spectrum checks. Case 429 has rates `(-.00048444,-2.82548,-93.11253)` and return period `10.64196533199` tau; case 1612 has `(.00010250,-3.45498,-46.69126)` and period `7.73800571098` tau. Independent source-form integrations close their periods within `4.09e-11`. The weak leading rates are therefore not promoted as chaos.

Two numerical problems are explicitly preserved. The first full-QR attempt lost contracting directions with one-unit QR spacing and produced no accepted spectra; the corrected run uses adaptive shorter intervals and a conditioning check. Separately, the initial focused-return driver stopped collecting after a source-algebra assertion, although workers independently saved 30 successful cases. The two rejected checks subtracted large reaction fluxes to obtain small net rates. A 70-digit Decimal source-field audit of all 32 initial states gives flux-scaled error below `1.92e-16` and complex-step Jacobian error below `1.26e-12`; raw residuals remain recorded. The two missing trajectories are completed in a separate run. The original failed driver record is not rewritten as a success.

Rho waves and complex oscillations already have substantial precedent. [Michaud et al. (2022)](https://doi.org/10.1083/jcb.202203017) addresses spatial Rho–actin behavior; [Tong, Xu and Wu (2023)](https://doi.org/10.1016/j.celrep.2023.112857) is a relevant coupled-signaling predecessor, whose complete equations were not recovered here. No claim of first Rho chaos follows from these searches.

## Two microbial populations with continuous cross-feeding

The source is [Ross et al., Nature Communications 16, 8919 (2025)](https://doi.org/10.1038/s41467-025-63986-y). The full five-state model is transcribed from [the authors' code](https://github.com/VenturelliLab/Ross_et_al_2025), pinned to commit `bde791b4c82cd3bebe8d2b34553d8e50c7b47d01`; author-file hashes are archived. States are two abundances `N1,N2`, two cross-fed resources `R1,R2`, and glucose `G`.

Each population's growth is the literal minimum of its resource-limited and glucose-limited Monod rates. Release depends on unused glucose-processing capacity. The calculation preserves these continuous, piecewise-smooth branches, with no smoothing, abundance clipping, periodic dilution or extra feedback. Time is `tau=D t`, with source time `t` in hours and positive dilution rate `D`.

There is an exact positive linear budget. In the author's parameter notation, let

```math
w=\left(q_{13}+\frac{q_{12}q_{23}}{q_{22}},\;
q_{23}+\frac{q_{21}q_{13}}{q_{11}},\;
\frac{q_{13}}{q_{11}},\;\frac{q_{23}}{q_{22}},\;1\right).
```

Then `B=w·(N1,N2,R1,R2,G)` satisfies `dB/dtau=B_in-B`. Positive initial data and positive parameters preserve the nonnegative orthant, so this budget bounds every component. The asymptotic budget hyperplane has four dimensions; the identity does not exclude chaos.

The source cycle control agrees between DOP853 and an independent stoichiometric Radau implementation within `2.39e-9` in scaled state over 40 tau. The source Jacobian check is below `2.85e-14`; the budget history residual is below `1.78e-14`.

A declared ten-parameter design uses 4,096 Sobol tuples and the source tuple, each with two swapped abundance starts: **8,194 lanes**, warmed for 200 and observed for 1,000 tau. The design varies dilution, three inlet concentrations, one maximum growth rate, glucose affinities, cross-feeding yields and paired resource affinities. Exact bounds and seeds are archived; this is not a fitted experimental ensemble.

The first implementation completes 5,552 lanes, with 1,766 step-underflow failures at limitation switches and 876 cases below the declared `1e-18` abundance threshold. Its normalized tangent error test over-resolves individual tangent components near zero. A separate rerun uses a vector-norm tangent tolerance while keeping state tolerances and the physical equations unchanged. It completes **7,318 lanes**, with **876 tiny-abundance cases unresolved and no step-underflow failures**. These thresholds do not prove extinction or exclude other dynamics.

The rerun retains **zero primary candidates** under the four-window `.005` rule. All 44 exploratory weak cases have final tangents aligned with the flow within `1.12e-15`. Independent 400-tau Radau checks of the source and seven strongest selected weak cases all give lag-one numerical returns. The source and three pilot IDs are checked first; four additional high-rate IDs are checked after the corrected rerun. A frozen protocol's inherited selection-description error is corrected in an explicit companion record; its actual job list and data are unchanged.

Cross-feeding chaos itself is already reported in [Shimada and Kaneko's 2026 preprint](https://arxiv.org/abs/2601.22564), with many populations and fast resource elimination. [McLaughlin et al.'s 2025 preprint](https://doi.org/10.21203/rs.3.rs-8021262/v1) treats two-strain chaos under serial feast–famine dilution. These differ from the continuous five-state Ross equations, but preclude a broad first-cross-feeding-chaos claim. Surplus-sharing feedback also has predecessors in [Ledder et al. (2020)](https://doi.org/10.1007/s12080-020-00464-9) and [Pfab et al. (2022)](https://doi.org/10.1093/conphys/coac026); no exact conjugacy is asserted.

## Motor–clutch dynamics: a model interpretation gate

[Shen and Zhang, Physical Review E 111, 014417 (2025)](https://doi.org/10.1103/PhysRevE.111.014417), [preprint](https://arxiv.org/abs/2402.18363), gives five equations for motor and clutch binding, displacements and substrate deformation. We implement the literal smooth Equation 6 / Appendix A1 and independently reproduce representative stable and periodic controls at `Nm=100,Nc=4,ks=kc`. The oscillatory control at `vu=1.5` agrees between DOP853 and Radau within `5.65e-9` over 4,000 units. A fourth control is stopped and remains unresolved.

The source's earlier piecewise force–velocity law and its later smooth Equation 6 are not interchangeable. Equation 6 continues the linear velocity branch beyond stall; roughly 44% of the sampled periodic control lies beyond stall. In that literal equation,

```math
\frac{d}{dt}(y+x_o)=v_u\left(1-\frac{k_my}{f_s}\right).
```

Any bounded nonconstant recurrent solution must therefore visit both sides of stall. Restricting it to the below-stall branch would change the result, not merely improve a numerical method. The substrate-equilibrium expression also omits a factor retained in the dynamical equation, and the printed force-scale unit is inconsistent with its defining scale. These differences are recorded before any broad chaos search. No claim is made against the source's experimental motivation, and no surrogate physical law has been silently substituted.

## Decision

The completed runs justify no new strange-attractor claim in these models. Their negative screens apply only to the declared numerical designs; unresolved cases remain unresolved. The [four-row motor geometry candidate](motor-four-geometry.md) is a separate result with different evidence. Further searches should start from the recorded source equations and prior-art gates, rather than reusing weak finite-time positives as discoveries.

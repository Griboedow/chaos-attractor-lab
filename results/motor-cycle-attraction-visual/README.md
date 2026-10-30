# Full-space approach to the five-cluster periodic orbit

This packet retains three sequential numerical illustrations around the certified eight-row periodic orbit. They share six deterministic generic directions and six reference phases; every declared start is retained under both integration settings. These are finite numerical experiments, not a basin-size calculation.

| Folder | Full physical perturbation norm | Outcome after six periods |
| --- | ---: | --- |
| This folder | 0.6 | None of the six passed the specified approach-to-reference gate. Five also failed the long pointwise solver comparison. No asymptotic destination is classified. |
| `local-0p001/` | 0.001 | All six approach the reference orbit and pass the integration comparison. Final full-state distances are approximately `3.89e-9–1.30e-8`. |
| `local-0p03/` | 0.03 | All six approach the reference orbit and pass the integration comparison. Final distances are approximately `1.32e-7–3.42e-7`. This is the displayed experiment. |

The first attempt used a larger offset to seek a visible transient; the two later amplitudes were chosen after inspecting the earlier outcomes. No individual trajectory was selected out of an ensemble for the final figure. The original driver version and its failed control outcome remain here. A completed integration is not the same as a passed solver-agreement gate: inspect the explicit flags in each `analysis.json`.

Each protocol records the initial states, parameter values, random seed, phase choices and source hashes before integration. The source row equations run in all 23 relative coordinates, without imposing pair equality. The physical perturbation has zero total displacement and Euclidean norm measured over all 24 `(c,s,X)` components. The exact parameter expressions of the certificate are converted to nearest floating-point values for these numerical integrations.

Each trajectory archive contains both complete solver histories, its initial state, parameters, sampled full-state distance to the reference orbit and the matched orbit phase. Distance uses a 12,001-point nearest-neighbor search followed by local dense-orbit phase minimization. The reference integration starts from the certified shooting-box midpoint; this midpoint and the numerical trajectory are not claimed to be an interval enclosure of the exact orbit.

The control uses DOP853 at maximum steps `0.04` and `0.02`, relative tolerances `2e-11` and `2e-12`, and absolute tolerances `2e-13` and `2e-14`. The displayed experiment's maximum difference over every component and stored time is below `8.62e-11`.

[Figure and interpretation](../../docs/attractor-geometry.md#approach-to-the-periodic-orbit) · [Attraction theorem](../../docs/motor-five-cycle-certificate.md) · [Final inventory](manifest.json).

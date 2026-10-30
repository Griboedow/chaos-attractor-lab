# Longer continuation of the eight-row motor candidate

This report extends the original [activity-path numerical evidence](motor-rows-activity-path.md) without selecting a new CPU trajectory. The full 23-dimensional CPU state, tangent basis, accumulated QR logarithms, divergence integral, flow logarithm and all recorded histories are restored from the frozen 4,000-observed-tau checkpoint. The additional segment uses the saved floating-point parameters, not newly rounded replacements.

## Full-dimensional CPU continuation

The input is [`cpu-N8-row69-extended.json`](../results/motor-rows/activity-path/cpu-N8-row69-extended.json) and its [NPZ checkpoint](../results/motor-rows/activity-path/cpu-N8-row69-extended.npz), at elapsed original tau `4199.999999998334`, including the original 200-tau warmup. The restored parameters are `(rho, r, b) = (2.934090909090909, 4.7551460477513405, 1.0)`. No new warmup or tangent reset is applied.

The continuation integrates the state and all 23 tangent directions using DOP853, relative tolerance `1e-11`, absolute tolerance `1e-13`, maximum integration step `0.03` and QR interval at most `0.2`. Every new 1,000-tau block saves its states, Poincare returns, endpoint basis, cumulative logarithms and separate block statistics. The [executed source](../research/continue_motor_activity_cpu.py), [terminal run record](../results/motor-activity-continuation/cpu/run-01/analysis.json), [independent analysis](../results/motor-activity-continuation/cpu/analysis/analysis.json) and [frozen CPU inventory](../results/motor-activity-continuation/cpu/manifest.json) retain the full audit trail.

The run completed all 16 new blocks: **20,000 observed tau**, with 200 additional warmup tau inherited from the original run. Wall time was 1,687.31 seconds. The six positive estimates remain separated from zero; the seventh estimate is consistent with the neutral flow direction at this finite duration.

| QR direction | Original 4,000 tau | Continued 20,000 tau |
|---|---:|---:|
| 1 | 0.213909 | 0.213134 |
| 2 | 0.185666 | 0.182374 |
| 3 | 0.157829 | 0.154626 |
| 4 | 0.114946 | 0.114598 |
| 5 | 0.080106 | 0.074331 |
| 6 | 0.028208 | 0.021640 |
| 7 | -0.000331 | -0.000058 |
| 8 | -0.021117 | -0.024626 |

The sixth estimate is positive in every new 1,000-tau window, with minimum `0.010174`. All 23 final estimates and all window values are in the analysis JSON. The final array retains 101,000 state samples and 58 complete restart checkpoints.

![CPU continuation and full spectrum](../results/motor-activity-continuation/cpu/analysis/cpu-spectrum.png)

### Checks independent of the leading exponents

The original histories are retained as byte-identical prefixes of the final arrays. The terminal run record and independent audit check all original samples, section crossings, checkpoints and old time windows. Initial and final endpoint checks compare the analytic Jacobian with complex-step differentiation, the normalized field with an independently solved source force constraint, and the analytic divergence with the Jacobian trace.

The transported flow direction is checked against the vector field after every QR interval. Its accumulated logarithm is also compared with the endpoint identity using the original flow norm at elapsed tau 200. This origin is read from the first CPU record: the extended record's `warm=0` describes its previous added segment and must not replace the original warmup.

The [fresh step control](../results/motor-activity-continuation/cpu/strict-control.json) starts from the exact inherited state and all 23 basis vectors. Over ten tau, the baseline reproduces the first 50 saved new states bitwise. Halving the maximum step to `0.015` and tightening the tolerances to `2e-13` and `2e-15` changes the states by at most `2.11e-13`, the final basis by `5.42e-13` and the logarithm sums by `4.71e-13`. This verifies the short numerical continuation; it does not give a global error bound for a chaotic trajectory over 20,000 tau.

A new scalar quadrature independently checks the [weighted-volume identity](motor-volume-contraction.md). With $\vartheta=7/8$, $a_8=-16+7(r-\rho)$ and $\Phi=\vartheta(r\sum_jc_j+\rho\sum_jX_j^2/2)$, the residual is

```math
R=\int\operatorname{div}f\,d\tau+\Phi_{\rm end}-\Phi_{\rm start}
-a_8\Delta\tau+\vartheta\int\sum_jV_j^2\,d\tau.
```

The speed-square quadrature begins at the restored checkpoint, so this check covers only the new 16,000-tau segment. The older segment has no saved speed integral. The QR determinant and integrated divergence are checked separately throughout.

The total new-segment weighted-volume residual is `4.07e-10`, against integrated divergence `-140570.45`; the largest local residual is `3.35e-14`. The accumulated QR/trace discrepancy per observed tau is `1.15e-13`. The independently transported flow estimate is `2.33072988149e-5`, agreeing with the endpoint telescoping value within `8.68e-16`.

Every one of the 28 physical row-pair distances is checked in each new block. The minimum pair RMS over all 16 blocks is `1.05613`, and in the final block it is `1.07421`; no pair is bitwise identical at the final state. Thus the saved samples do not show the pair-synchronization capture found in some separate GPU lanes. As an additional finite diagnostic, full labeled-row distances between the last 1,000 Poincare returns are checked at lags 1 through 64. The smallest RMS is `2.47031` at lag 2, with smallest individual distance there `0.59925`; this does not establish aperiodicity or exclude longer returns.

### Provenance and reproduction

The first preflight stopped before integration because the original short-run JSON had been normalized from CRLF to LF since its historical hash was recorded. Both current bytes and the exactly reconstructed historical CRLF bytes are preserved; the latter reproduce the recorded hash. The failed [preflight record](../results/motor-activity-continuation/cpu/preflight-01.json) and its executed source remain in the evidence inventory. No old input file was rewritten.

Run these commands from the repository root with NumPy, SciPy and Matplotlib installed. Use fresh output paths. The distant cutoff below avoids reusing the historical run's absolute wall-clock deadline; it does not change the integration settings.

```console
python research/continue_motor_activity_cpu.py --output outputs/motor-cpu-recheck --cutoff-utc 2099-01-01T00:00:00+00:00
python results/motor-activity-continuation/cpu/strict-control.source.py --run outputs/motor-cpu-recheck --output outputs/motor-cpu-strict-control.json
python results/motor-activity-continuation/cpu/analyze.source.py --run outputs/motor-cpu-recheck --output outputs/motor-cpu-analysis
```

All reported exponents are finite-time numerical estimates in original tau. These checks strengthen the retained candidate evidence, but do not prove asymptotic hyperchaos, an attracting invariant set, a basin size or scientific novelty. The formal global attractor and zero-volume results are separate [analytic statements](motor-global-geometry.md).

## The unchanged 32-member GPU cohort

Both GPU continuations start from the same saved states and carried tangents at absolute `tau=7100` in [N8-lifetime-history.npz](../results/motor-rows/activity-path/N8-lifetime-history.npz). They retain all 32 source IDs, including the previously subthreshold lane49. There is no new warmup, pruning or reselection. Each setting adds forty 500-tau blocks, each containing four 125-tau windows, reaching absolute `tau=27100`.

The [main run](../results/motor-activity-continuation/gpu-main/summary.json) uses relative tolerance `1e-10` and step factor `0.5`; the [strict run](../results/motor-activity-continuation/gpu-strict/summary.json) uses `2e-11` and `0.25`. The original scanner normalizes the length of the carried diagnostic tangent at each block, while preserving its direction. All 2,560 lane-block integrations complete. The same **31 of 32** members pass every new block mean threshold `0.005` and every window threshold `0.001` under both settings. This is a result for a selected historical cohort, not a basin fraction.

The [independent GPU audit](../results/motor-activity-continuation/gpu-audit-recheck/analysis.json) passes 268 checks of source hashes, complete durations, exact state/tangent carry, unchanged IDs and parameters, and final checkpoints. It also records physical row distances at every saved block endpoint. The earlier audit is preserved; the recheck adds these endpoint histories.

| Source lane | rho | Main mean over new 20,000 | Strict mean over new 20,000 | Smallest main/strict 125-tau window |
|---|---:|---:|---:|---:|
|36|2.0090909091|0.269451|0.283155|0.081558 /0.095859|
|42|2.1772727273|0.247773|0.249327|0.098394 /0.113405|
|46|2.2613636364|0.043557|0.043289|0.006027 /0.012945|
|49|2.3454545455|0.00000114|0.00000114|−0.003054 /−0.003054|
|60|2.6818181818|0.222932|0.255197|0.097490 /0.133793|
|69|2.9340909091|0.210773|0.212056|0.174320 /0.169146|

![All retained GPU histories and endpoint pair structure](../results/motor-activity-continuation/gpu-audit-recheck/cohort-continuation.png)

Positive tangent growth requires an additional qualification here. Ten positive members already have synchronized pairs at the initial checkpoint. The strict lane60 also first has a saved row distance below `1e-10` at absolute `tau=9100` and ends with three bitwise-equal pairs; the main lane60 has no such near pairs at its saved endpoints. Both remain positive. Therefore the two settings reproduce the overall persistence count without proving identical long-time cluster behavior. The representative lane69 has no near or exact pair at any saved endpoint in either setting, and its separate CPU history supplies the denser all-pair check above.

The [continuation driver](../research/continue_motor_activity_gpu.py) and [analysis program](../research/report_motor_activity_gpu.py) refuse to overwrite outputs. Reproduce the two integrations from the repository root with fresh folders:

```console
python research/continue_motor_activity_gpu.py --setting main --output outputs/activity-gpu-main
python research/continue_motor_activity_gpu.py --setting strict --output outputs/activity-gpu-strict
```

## The subthreshold member is a periodic saddle

Lane49 illustrates why neither a near-zero carried exponent nor a fresh positive exponent classifies an attractor. At its archived checkpoint, three row pairs are bitwise equal and the fourth differs by about roundoff. A separate full-state sample and weighted eleven-dimensional shooting calculation locate a periodic state with pairs `(0,3)`, `(1,6)`, `(2,5)` and `(4,7)` and period `15.65798354606`.

The [periodic calculation](../results/motor-activity-continuation/lane49/periodic/analysis.json) has closure residual `2.08e-13`. Independent DOP853 and Radau quotient monodromy matrices differ by at most `8.59e-13`. All ten nonflow quotient multipliers have modulus below `0.748`. However, the full 23-dimensional monodromy has **three unstable, one flow and nineteen stable directions**. This is a numerically resolved periodic saddle, not an interval certificate or an attracting periodic orbit.

An [independent contrast-block control](../results/motor-activity-continuation/lane49/transverse-control.json) finds one multiplier approximately `−1.72783956727` in each of the first three pair blocks; their expanding exponents are about `0.03492607`. The fourth pair block is stable. DOP853 and Radau block matrices agree within `1.48e-13`. A one-third-period shift cyclically permutes the first three pairs while leaving the fourth fixed, with state residual `2.82e-13`. The retained [provenance record](../results/motor-activity-continuation/lane49/transverse-control-provenance.json) identifies the statements executed through standard input and the source hashes recorded afterward.

The [same-state restart control](../results/motor-activity-continuation/lane49/restart-control/analysis.json) keeps the exact archived physical state and changes only the diagnostic tangent. Over 500 tau, the carried tangent gives mean `−0.00142063`, whereas a fresh generic tangent gives `0.03387877`. The latter detects the periodic saddle's hidden transverse instability; it is not evidence of chaotic motion of that unchanged state. Two separate physical initial perturbations of pair `(0,3)`, of sizes `1e-7` and `1e-9` in its c half-difference, reach a full pair-distance threshold `1e-3` after about 214.07 and 339.53 tau. These are declared new trajectories, not recovery of digits already lost at the saved checkpoint.

The original 31/32 screening count remains unchanged. The new result explains the subthreshold member's numerical capture and shows why cluster-restricted dynamics need independent transverse checks.

## Conditional growth and physical perturbations

Three preselected clustered endpoints receive a separate weighted-cluster calculation. The restrictions retain their true fractions of the eight rows. The measured initial projection errors are zero for lane36, `1.70e-14` for lane42 and `5.59e-15` for lane46. These explicitly projected diagnostic trajectories do not replace any original cohort member.

The [conditional driver](../research/check_motor_activity_clusters.py) evolves one quotient tangent and a complete 3-by-3 tangent basis for each paired-row contrast, with 200 warmup tau and 2,000 observed tau. Independent checks embed each contrast block into the full 23-dimensional Jacobian; the maximum discrepancy across these cases is `1.78e-15`. Every 500-tau window and all three rates per block are retained.

| Source lane | Cluster sizes | Quotient dimension | Leading quotient rate | Leading conditional rates of the paired clusters |
|---|---|---:|---:|---|
|36|2,2,1,2,1|14|0.276729|0.020929,0.000161,−0.060727|
|42|2,2,2,2|11|0.269667|0.009817,0.012681,−0.004893,0.000517|
|46|1,2,2,2,1|14|0.044531|−0.092096,−0.090095,−0.077198|

For [lane36](../results/motor-activity-continuation/cluster-controls/lane-36/analysis.json) and [lane42](../results/motor-activity-continuation/cluster-controls/lane-42/analysis.json), the window estimates change sign substantially. Their small averaged conditional rates do not settle asymptotic transverse stability. For [lane46](../results/motor-activity-continuation/cluster-controls/lane-46/analysis.json), every leading conditional rate is negative in every one of the four windows. This is a more promising finite diagnostic of transverse attraction, while the quotient remains irregular with positive growth.

A further [declared perturbation protocol](../results/motor-activity-continuation/full-space-perturbations/protocol.json) tests the original physical endpoints 36, 42, 46 and 49. For each, it retains the unchanged state with a fresh diagnostic tangent and four new states with physical perturbation norms `1e-7` and `1e-9`, in both signs along one fixed generic direction. The direction has zero total X displacement and unit norm in all 24 physical coordinates. The reduced 23-state representation and the actual post-rounding norms are saved. Seed 2026103049 fixes the four directions; this small targeted test is not random basin sampling.

All 20 cases complete 1,000 warmup tau plus 4,000 observed tau at the strict GPU settings. The [independent audit](../results/motor-activity-continuation/perturbation-audit/analysis.json) passes 89 checks of hashes, exact starts, physical norms, parameters and full durations.

![Every declared physical perturbation outcome](../results/motor-activity-continuation/perturbation-audit/perturbation-outcomes.png)

All four perturbed lane46 starts keep every observed window positive, with overall means between 0.04035 and 0.04793. They end near the original three paired groups. Lane36 and lane42 perturbations also remain positive, but some lane42 outcomes have different pair assignments. These observations describe finite endpoint behavior and cannot establish basin size or asymptotic normal attraction.

For lane49, the unchanged physical state with a fresh tangent retains a positive mean near 0.03515, consistently with its known periodic-saddle instability. All four physically perturbed starts fail the original screening rule by the last 1,000-tau window; their means over the entire 4,000-unit observation remain positive because of the earlier transitions. Their final windows lie between approximately −0.00058 and 0.00035. Separate return and Floquet calculations are required to classify their destinations. Positive finite-time means alone would miss this change.

Reproduce the separate diagnostics with unused output folders:

```console
python research/check_motor_activity_clusters.py --output outputs/activity-cluster-controls --lanes 36 42 46
python research/check_motor_activity_perturbations.py --output outputs/activity-physical-probes
python research/report_motor_activity_perturbations.py --input outputs/activity-physical-probes --output outputs/activity-probe-audit
```

The [completed-packet inventory](../results/motor-activity-continuation/publication-02/manifest.json) passes 119 additional provenance and completion checks and hashes 299 archived files. It preserves the publication text and reproduction sources separately from any later additions. The first inventory attempt failed on relative versus absolute path handling after the scientific checks passed; its source and failure record remain in `publication-01`. Local Markdown paths, math fences and the three PNG figures were checked; browser rendering was unavailable in this session.

## Completed follow-ups: a cluster candidate and a periodic destination

A [complete conditional spectrum for lane46](motor-lane46-spectrum.md) now resolves all 14 quotient directions and all 9 pair-splitting directions over 4,000 observed tau. It has one positive quotient estimate near 0.0450, one near-neutral flow estimate, and 21 negative estimates. Every contrast leader is negative in every 500-tau window. A separate tighter 1,000-tau return diagnostic finds no close periodic closure at the tested section lags 1 through 64. These finite results support a cluster-chaos candidate with transverse contraction; they remain short of an asymptotic attraction theorem.

All four near-zero lane49 probe outcomes also receive a [separate destination study](../results/motor-activity-continuation/probe-destinations/README.md), including 500 additional tau in the unprojected 23-dimensional field, section-return tests, shooting and full monodromy matrices. Three shoot back to symmetry-related copies of the original periodic saddle, matching after phase and row permutation to below `1.7e-13` RMS. Thus their near-zero final windows did not establish stable capture.

The `+1e-7` case instead reaches a **distinct five-cluster periodic candidate**, with groups `[[0,3],[1,2],[4],[5,6],[7]]` and period approximately 25.00168437590. Its full numerical spectrum has one flow multiplier and 22 stable multipliers; the largest nonflow modulus is about 0.12748229. A completely independent 23-dimensional Radau state calculation and all 23 variational columns confirm the spectrum: matrix disagreement is below `1.73e-11` and matched eigenvalue disagreement below `7.41e-13`. This is strong numerical evidence of a periodic attractor at the same parameter as the original saddle, with no interval or novelty claim in this numerical packet.

The [destination manifest](../results/motor-activity-continuation/probe-destinations/manifest.json) freezes all four outcomes, sources, refined phase comparisons and numerical symmetry tests. Any later interval theorem is a separate proof package. The initial `publication-02` inventory remains unchanged and deliberately excludes these later follow-ups.

# Overrelaxed Sinkhorn numerical packet

This packet retains all 576 declared inputs, every positive pilot flag and
both continuations of all 40 flags to iteration 5,530,000. It also contains
four full finite-window spectra, four arbitrary-precision runs, 24 state
controls and implementation checks. No new attractor portrait is promoted.

[Report, equations and limits](../../docs/sinkhorn-screen.md) ·
[Research index](../../README.md#research-index)

## Retained evidence

| Location | Contents |
| --- | --- |
| [Pilot](pilot) | All matrices, marginals through their defining kernel, starts, relaxation values and outcomes; 192 cases for each size |
| [First continuation](followup) | Twenty 25,000-step blocks per flag, including every block endpoint |
| [Long continuation](extension) | Twenty 250,000-step blocks per flag, including every block endpoint |
| [Spectra](analysis/spectra.json) | Four complete QR spectra and blockwise volume-sum errors |
| [Recurrence](analysis/recurrence.json) | All 40 final-window recurrence checks, approximate cycle multipliers and independent 100-digit return residuals |
| [Controls](analysis/controls.json) | Ordinary-update, near-origin and perturbed-endpoint controls; ordinary-update growth rates are not assessed |
| [Case 55 precision](precision01) | 70- and 100-digit positive finite windows before later periodic capture |
| [Case 64 precision](precision02) | 70- and 100-digit finite windows from the long-run checkpoint |
| [Checks](checks) | Independent probability-domain update, finite differences, determinant and inverse checks |
| [Code](code) | Executed source snapshots and packaging scripts |
| [External arrays](arrays.json) | Eleven optional state files, 46,513,332 bytes, with hashes, shapes and exact coverage |
| [Original summaries](original-summaries.json) | Original JSON hashes and value-preserving compacting or splitting into retained files |
| [Sources](sources.json) | Retrieval dates and hashes of external primary-source HTML copies |
| [Workflow](workflow.json) | Completed stages, development failures and their corrections |
| [Manifest](manifest.json) | Retained byte hashes and dependencies |

Case identifiers are local to a matrix size. The final-window files are
explicitly partial histories. The pilot retains endpoints and at most its
last 512 states; each continuation retains all states of its last 25,000
steps. Full earlier histories were not saved. Arbitrary-precision observed
states are stored as binary64 exports for summaries; the final high-precision
state is retained as decimal strings. Nothing is silently thinned under a
full-trajectory filename.

## Verification

From the repository root:

```bash
python research/verify_sinkhorn_screen.py
python research/verify_sinkhorn_screen.py --numerical --workdir D:/Codex/sinkhorn-recheck
```

The first command checks retained hashes, case accounting, lineage, reported
spectra and recurrence counts. The optional numerical check replays the
independent small-model controls in a new external directory. It does not
repeat the long campaign or prove asymptotic chaos.

## Reproduction

The recorded environment is Python 3.12.14, NumPy 2.5.4, SciPy 1.18.1,
Numba 0.68.1, mpmath 1.3.0 and Matplotlib 3.11.2. See
[environment.json](environment.json). Copy the files in `code/` into a new
external run directory, change to that directory and keep numerical outputs
there. Set `NUMBA_CACHE_DIR`, `MPLCONFIGDIR`, `TEMP` and `TMP` to external
locations. The survey and continuation drivers use four Numba threads;
use one BLAS thread to avoid nested parallelism.

Run the commands sequentially; later stages require completed parent files:

```bash
python check_model.py
python survey.py --output pilot02
python followup.py --pilot pilot02 --output followup01
python extend.py --input followup01 --output extension01
python spectrum.py --input extension01 --output spectra02
python analyze.py --input extension01 --output recurrence02.json
python precision.py --input followup01/n4.json --output precision01 --case 55
python precision.py --input extension01/n4.json --output precision02 --case 64
python controls.py
python check_actual.py
python diagnostic.py --output growth-diagnostics.png
```

The initial pilot source uses deterministic seeds, but chaotic long histories
need not agree pointwise across numerical libraries and platforms. Compare
declared inputs, exact structural identities, completed-window statistics and
capture classifications. A different capture time is scientifically relevant
and should be retained, not replaced by the archived trajectory.

The [report](../../docs/sinkhorn-screen.md#prior-art-and-equivalence-scope)
distinguishes known nonconvergence, finite growth, attracting dynamics and
historical novelty. None of these terms substitutes for the next.

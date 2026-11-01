# Hydrogel refinement and variational checks

[Results and interpretation](../../docs/hydrogel-variational-checks.md) ·
[Earlier thinner-rod screen](../hydrogel-slender-screen/README.md).

Twenty declared runs test endpoint restarts, 64/128-segment refinement,
time steps 0.006/0.003, three tangent directions, a rotating control and
longer continuation. Nineteen complete; the field-130 continuation fails
after developing subdiameter clearances. Positive finite growth at field
124 accompanies substantial slow drift. No chaotic attractor or historical
priority is established.

## Contents

| Record | Contents |
| --- | --- |
| `runs/*.json`, `runs.json` | Every integration report, including pilot and failure |
| `workflow.json`, `code/replay.py` | Declared seed graph and executable run sequence |
| `code/` | Executed solver, restart, tangent and diagnostic sources; old assembly retained |
| `checks/` | Derivative, symmetry, BDF-history, remeshing and assembly controls |
| `analysis/` | Full retained compact recurrence scans, clearance counts, drift windows and bounded source search |
| `arrays.json` | Hashes, shapes and external paths of complete trajectories and checkpoints; original seed identities |
| `environment.json`, `manifest.json` | Recorded environment and retained-file/dependency hashes |

The long-lag recurrence scan retains every tested lag and the refined local
minima. It is not a period certificate. Clearance records state the temporal
sampling and material-neighbor cutoff. Restart reports distinguish the last
saved sample from the last accepted step. All full numerical arrays remain
outside Git; no thinned trajectory is substituted under an existing name.

## Reproduction

Use Python 3.12 with NumPy, SciPy, Numba and Matplotlib; exact observed
versions are in `environment.json`. Set single-threaded BLAS and external
temporary, bytecode and Numba cache directories. From the repository root:

```sh
python research/verify_hydrogel_variational.py
python research/verify_hydrogel_variational.py --numerical --workdir D:/Codex/hydrogel-algebra-replay
```

The numerical verifier copies frozen code to a new external directory and
repeats synthetic assembly/remeshing, complete-step tangent, symmetry,
two-level QR and analytic relative-period checks. It does not require the
large trajectories, replay long integrations or certify chaos.

For a complete integration replay, first recreate field 117 through time 96
from the [slender packet](../hydrogel-slender-screen/README.md#reproduction)
and the 64-segment rotating control through time 100 from the
[spatial packet](../hydrogel-spatial-controls/README.md#reproduction).
Copy all files in this packet's `code/` to a new external code directory,
then run:

```sh
python D:/Codex/hydrogel-replay-code/replay.py --base117 D:/Codex/slender-replay/lambda127-field117 --rotation D:/Codex/spatial-replay/rotation-n64 --output D:/Codex/hydrogel-full-replay
```

The destination must not exist. This runs all 20 declared cases sequentially;
the anticipated failed field-130 continuation retains its accepted states.
Numerical failure times can vary by platform. The script stops if another
run unexpectedly fails. Recompute diagnostics by passing the desired run
folders to `clearance.py`, `recurrence.py` and `stationarity.py` with explicit
external output paths. For example:

```sh
python D:/Codex/hydrogel-replay-code/clearance.py D:/Codex/hydrogel-full-replay/spectrum124-extension --whole-history --output D:/Codex/hydrogel-full-replay/clearance-extension.json
python D:/Codex/hydrogel-replay-code/recurrence.py D:/Codex/hydrogel-full-replay/spectrum117-n64 --upper 20 --samples 991 --output D:/Codex/hydrogel-full-replay/long-returns.json
```

Original campaign checks that name local seed files are retained verbatim.
Use the portable verifier above for synthetic checks; those historical
scripts require their declared seed data. `compare_assembly.py` and plotting
scripts expect run folders beside their copied source. Keep those outputs
outside the checkout. The diagnostic figure is generated from the retained
reports and drift summaries by `render_diagnostics.py`; it is not an attractor
portrait.

[Research index](../../README.md#research-index) · [Storage policy](../../docs/data-storage.md)

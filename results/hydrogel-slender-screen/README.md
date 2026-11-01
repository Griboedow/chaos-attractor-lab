# Slender hydrogel filament screen

[Model, symmetry reduction and interpretation](../../docs/hydrogel-spatial-controls.md#circular-section-symmetry-and-thinner-filaments) ·
[Earlier spatial controls](../hydrogel-spatial-controls/README.md).

Three declared runs use slenderness 127, circular cross-section, fields
117, 130 and 150, 64 segments and BDF2 step 0.006. Two complete the requested
time 96; field 150 fails Newton convergence at 52.620. Both stronger fields
have sampled nonlocal clearances below the nominal diameter. The field-117
control is still unresolved in the tested late-window returns. No new
chaotic attractor or historical-priority claim follows.

## Records

| Record | Contents |
| --- | --- |
| `runs/*.json` | All three integration reports, including failure |
| `code/symmetry.py`, `symmetry-check.json` | Circular-section rotation action, global angular comparison and numerical controls |
| `code/returns.py`, `returns.json` | Full-state section returns with and without axial rotation |
| `rotation-control.json`, `control-input.json` | The same comparison on the previously recorded rotating control; input hashes |
| `code/clearance.py`, `clearance.json` | Finite sampled segment-distance screen and its explicit cutoff |
| `code/preview.py` | External diagnostic plot; not an attractor portrait |
| `arrays.json` | Hashes, shapes and locations of the three unthinned external histories |
| `source-audit.json`, `manifest.json` | Bounded source search and retained-file/dependency hashes |

The solver is unchanged from the earlier packet. Its hash is recorded in
each run and in the manifest. The minimum-over-lags diagnostic has different
pair counts at different lags; all counts and values are retained. It is
not a statistical test or a Lyapunov estimate. The rotating control has only
one height-return pair in its late window. Reported failures retain accepted
saved samples only; no failed Newton trial is appended as a trajectory state.

## Reproduction

Use Python with NumPy, SciPy and Numba, single-threaded BLAS, and external
temporary/cache directories. The recorded environment is Python 3.12.14,
NumPy 2.5.4, SciPy 1.18.1 and Numba 0.68.1. Run from the repository root:

```sh
python research/verify_hydrogel_slender.py --numerical
python research/hydrogel_spatial.py --segments 64 --dt 0.006 --duration 96 --slenderness 127 --eta 1 --chi 117 --output D:/Codex/slender-replay/lambda127-field117
```

Repeat the simulation command with fields 130 and 150 and distinct output
directories. The solver rejects existing destinations. To repeat the
postprocessing, copy the four files in `code/` to `D:/Codex/slender-replay/`
and copy `research/hydrogel_spatial.py` there as `spatial_fast.py`. Then run
`symmetry.py`, `returns.py`, `clearance.py` and optionally `preview.py` in
that external directory. These scripts preserve the full arrays; their
compact JSON outputs can be compared with this packet.

For the earlier rotating control, reproduce its 64-segment run from the
[spatial packet](../hydrogel-spatial-controls/README.md#reproduction), then
pass its folder explicitly to `returns.py` with
`--output D:/Codex/slender-replay/rotation-control.json`.
The verifier checks retained hashes, counts and numerical consistency;
`--numerical` replays the symmetry controls without writing into the checkout.
It does not replay long integrations or certify continuum dynamics.

Generated arrays and preview images remain outside Git. No new portrait,
PDF or reduced replacement trajectory is included.

[Research index](../../README.md#research-index) · [Storage policy](../../docs/data-storage.md)

# Hydrogel planar controls

Independent numerical controls for the planar restrictions of
[Boiardi et al. (2026)](https://arxiv.org/html/2604.09486v1).
[Report and equations](../../docs/hydrogel-planar-controls.md).

The packet contains two spatial discretizations of the linearized problem,
four nonlinear periodic controls, all eight exploratory cases and three longer
follow-ups. Six scout trajectories have sampled centerline self-intersections;
the two others close to numerical periodic orbits in longer checks. No new
chaotic attractor or three-dimensional stability result is asserted.

## Contents

| Record | Scope |
| --- | --- |
| `linear.json` | Hermite and angle thresholds at four mesh sizes; passive, drag and linearization controls |
| `runs/*.json` | Complete compact records of four controls, eight scouts and three extensions |
| `inspection.json` | Sampled proper self-intersections for all fifteen runs |
| `arrays.json` | Sizes, shapes and SHA-256 hashes of external, unthinned numerical histories |
| `code/` | Exact executed model, eight-case driver and crossing inspection source |
| `source-audit.json` | Source version and bounded literature-search scope |
| `manifest.json` | Hashes of the retained packet and reusable model |

Arrays and the source HTML are stored outside Git on D:. No PDF, portrait,
thinned replacement trajectory or browser screenshot is included in this
numerical packet. Initial states, controls, tolerances, time intervals and
solver versions remain recorded. Root-refined event returns are distinct from
the uniformly sampled histories used for crossing detection.

## Reproduction

The recorded environment was Python 3.12.14, NumPy 2.5.4, SciPy 1.18.1 and
Numba 0.68.1. Install `numpy scipy numba` in an isolated environment. Set
`NUMBA_CACHE_DIR`, `TEMP` and `TMP` to an external working directory. Setting
`OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1` avoids nested linear-algebra
threads when running several cases.

From the repository root, use a new external destination for each command:

```sh
python research/verify_hydrogel_planar.py --numerical
python research/hydrogel_planar.py linear --output D:/Codex/hydrogel-replay/linear
python research/hydrogel_planar.py run --segments 32 --duration 60 --chi 100 --rtol 1e-8 --method Radau --output D:/Codex/hydrogel-replay/control
python research/hydrogel_planar.py run --segments 64 --duration 600 --chi -33 --rtol 1e-8 --method Radau --output D:/Codex/hydrogel-replay/negative
python research/hydrogel_planar.py run --segments 64 --duration 400 --slenderness 127 --eta 1 --chi 117 --rtol 1e-8 --output D:/Codex/hydrogel-replay/slender
```

The verifier checks retained hashes and summary consistency. `--numerical`
also reruns both linear discretizations and the independent drag quadrature;
it does not replay all nonlinear histories or certify continuum dynamics.
The standalone model works from another current directory when invoked by
absolute path. Output directories must not already exist.

To repeat the eight-case screen, copy the three files in `code/` into an empty
external directory and run `scout.py` there. This original driver writes
`scout-01/` beside its own source and uses three worker processes. The model
accepts the individual parameters listed in every retained run, so the same
cases can also be run serially with the command-line interface. The frozen
`inspect_runs.py` expects the original run-folder names beside its source;
it inspects the available completed runs and does not reconstruct missing data.

The source article's full spatial computation uses additional director and
twist equations. Neither these scripts nor the planar periodic controls are a
replacement for that three-dimensional validation.

[Research index](../../README.md#research-index) ·
[Data storage](../../docs/data-storage.md)

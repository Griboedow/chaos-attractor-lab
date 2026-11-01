# Hydrogel spatial controls

[Report, equations and mesh artifact](../../docs/hydrogel-spatial-controls.md) ·
[Earlier planar validation](../hydrogel-planar-controls/README.md).

This packet records a spatial segment approximation to the published
morphoelastic rod, including its algebraic twist balance. It contains all
twelve declared runs: ten complete and two fail Newton convergence. Numerical
periodic beating and hook-shaped rotation are controls; a coarse bent
equilibrium is a discretization artifact. The completed stronger-field scout
has sampled clearances below the nominal filament diameter. No new chaotic
attractor, physical prediction or priority claim is made.

## Retained records

| Record | Scope |
| --- | --- |
| `code/spatial.py` | Original complete-difference implementation |
| `code/spatial_fast.py` | Executed faster implementation; same equations |
| `code/*check*.py`, `code/*compare.py` | Executed algebra, Jacobian, equilibrium and temporal checks |
| `algebra.json`, `fast-check.json` | Independent energy/drag checks and comparison of Newton implementations |
| `planar-comparison.json` | Full histories and periods compared with the earlier adaptive planar solver |
| `spurious-equilibrium.json` | Explicit coarse bent equilibrium and comparison with the observed state |
| `runs/*.json` | Parameters, integration completion/failure, residuals and finite summaries |
| `analysis.json`, `rotation-summary.json` | Sampled section returns and late angular fits |
| `clearance.json` | Sampled nonlocal polygonal distances and explicit test limits |
| `arrays.json` | External, unthinned trajectory hashes and shapes |
| `source-audit.json`, `manifest.json` | Source scope and retained-file hashes |

The raw run field `max_hinge_angle` includes the final Newton trial when that
trial fails. A failed run's saved history contains accepted samples only.
`completed_time` is the last saved time. Section returns use cubic
interpolation of saved quaternion and curvature histories, not interval
enclosures or the integrator's exact event interface. Quaternion signs are
identified when comparing orientations. Angular-fit periods for the rotating
control are distinct from root-refined full-state returns.

Large arrays and diagnostic previews remain on D:. No PDF, new attractor
portrait or thinned replacement array is committed.

## Reproduction

The calculations used Python 3.12.14, NumPy 2.5.4, SciPy 1.18.1 and Numba
0.68.1. Install `numpy scipy numba` in an isolated environment and keep Numba
cache and temporary directories outside the checkout. Single-threaded BLAS
was used. Every destination below must be new:

```sh
python research/verify_hydrogel_spatial.py
python research/hydrogel_spatial.py --segments 32 --dt 0.002 --duration 60 --output D:/Codex/hydrogel-spatial-replay/beating
python research/hydrogel_spatial.py --segments 64 --dt 0.004 --duration 100 --eta 1 --chi -33 --output D:/Codex/hydrogel-spatial-replay/rotation
python research/hydrogel_spatial.py --segments 64 --dt 0.004 --duration 60 --chi 150 --output D:/Codex/hydrogel-spatial-replay/strong-field
```

The model can be invoked by absolute path from another current directory.
The reusable source is byte-identical to `code/spatial_fast.py`.
The verifier checks retained hashes, completion counts and consistency of the
compact evidence. Its optional `--numerical` repeats an energy derivative
check, both BDF Jacobian checks and the constructed discrete equilibrium;
it does not replay twelve integrations or establish continuum convergence.

For the original implementation, copy `code/spatial.py` to an external working
directory and use the same command-line arguments. The faster implementation
is appropriate for new runs; both implementations solve the same mechanical
residual. The frozen exploratory postprocessors use original run-folder names
beside their source. `check_algebra.py` and `planar_compare.py` also refer to the
recorded earlier D: checkout/cache; these historical entry points are preserved
as executed. The portable verifier supplies current standalone algebra checks.

Original trajectories are optional replay outputs, not hidden simulation
inputs. Keep failed and completed runs distinct: a nonlinear solver failure
does not show that the source model escapes or lacks an attractor.

[Research index](../../README.md#research-index) ·
[Storage policy](../../docs/data-storage.md)

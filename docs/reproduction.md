# Reproduction

Use the instructions and environment recorded for the selected study. The repository contains several numerical pipelines; one environment is not a validation of every result. Write new outputs to a separate directory.

Large computed arrays are optional cache data. See [data storage](data-storage.md)
for a shallow clone, selected-payload recovery and the `--materialize` option
required by older verification scripts. Sources, small inputs, compact summaries
and certificates remain tracked.

## Selected verification entry points

| Study | Instructions |
| --- | --- |
| Icosahedral flow | [Search and absorbing bounds](../research/icosahedral_search/README.md), [equilibria](icosahedral-structure.md), [numerical returns](../research/icosahedral_returns/README.md) |
| RAAR symbolic dynamics | [Certificate and manuscript](../papers/raar-symbolic-dynamics/README.md) |
| Scalar acceleration memory | [Reproducible search](../research/unexpected_search/README.md) |
| Adaptive solvers | [RK45](../papers/rk45-controller/README.md), [RK23](../papers/rk23-controller/README.md) |
| Optimizer memory | [Yogi](../papers/yogi-memory/README.md), [Adamax](adamax-geometry.md) |
| Motor periodic attractor | [Five-cluster certificate](motor-five-cycle-certificate.md), [proof implementation review](motor-proof-review.md) |
| Presentation figures | [Rendering and date conventions](../tools/README-portraits.md) |

## Recorded environments and runtime checks


The later motor modules were also checked in the updated independent public clone at `ad83bf2`: N=3/N=8 GPU versus CPU controls passed under isolated Python, with no Chua imports or changes to archived results. [Additional packaging record](../provenance/motor-public-clone-validation.json).

The published five-cluster proof tree at `bcd5119` was separately exported from the fetched public branch. All 24 manifest hashes matched; its algebra controls and exact symmetry enumeration passed in the separate proof environment. Local links were checked, while browser rendering was unavailable for these latest reports. [Publication verification](../provenance/motor-session-publication-20261030.json).

**Public-clone check:** a fresh clone of this repository at `a514509`, with a newly installed environment, passed four GPU/CPU runtime controls, the optional DCR control and `pip check`. The checks ran outside both repositories with isolated Python, imported no Chua module and left the archived evidence unchanged. [Validation record](../provenance/public-clone-validation.json), [file-copy audit](../provenance/publication-copy-integrity.json). This verifies the tested packaging paths, not every scientific result or platform.

From the project root, create an environment with `python -m venv .venv`. Activate it using `source .venv/bin/activate` on Unix-like systems or `.venv\Scripts\Activate.ps1` in PowerShell, then run `python -m pip install -r requirements.txt`. Python 3.14.3 was used for the recorded work. A separate fresh environment passed installation, `pip check` and four extracted GPU/CPU controls at source snapshot `8a3564c`; those tested modules and the runtime extraction are unchanged in this refreshed copy. The [validation record](../provenance/packaging-validation.json) retains those checks, and [requirements-lock.txt](../requirements-lock.txt) records the observed environment. This does not certify every model, every platform or a future public clone. GPU scripts require a working OpenCL GPU driver with FP64 support. They do not silently substitute CPU execution.

The optional archived DCR source-solver controls additionally require `python -m pip install -r requirements-dcr.txt` after installing the main requirements. This adds pinned `scikit-sundae`/SUNDIALS bindings. Version 1.2.0 passed a separate fresh-environment dependency check and an isolated CPU control of duration 0.01. Its maximum state discrepancy from the archived fine control was 1.16e-10; [optional-dependency provenance](../provenance/optional-dcr-packaging.json) retains the original 3cf58d0 snapshot and confirms the relevant files are unchanged here. This is a packaging check, not a new scientific result. The observed lock includes this optional package. It was not part of the four-module base-runtime validation and no DCR mass survey is supported. See [the DCR reproduction instructions](../docs/dcr-research.md) before starting its bounded controls.

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

Plotting commands write figures. Keep the archived evidence in a clean checkout and use a disposable copy when reproducing figures or simulations. The [results index](../results/README.md) and individual reports identify settings, seeds, traces and exact model versions. Many simulation scripts overwrite named result files; read their arguments before running.

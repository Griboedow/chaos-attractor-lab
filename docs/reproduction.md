# Reproduction

Run commands from the repository root. Each study records its own equations, parameters, initial states and software versions. Use an external output directory or a disposable checkout for scripts that write to fixed paths.

## Check the two-ring attractor

Case 13745 has two annuli that map strictly into each other. To repeat the interval proof, use a Python 3.11 virtual environment and the certificate committed with the study. The verifier runs on a CPU and reads the exact coefficients from the repository.

```sh
python -m pip install python-flint==0.9.0
python research/wow_search_20261122/symmetry/certify.py --verify --case 13745
```

The verifier checks 392 boxes in the inner annulus and 132 in the outer annulus, then compares the bounds with the stored certificate. [Proof and full portrait replay](../results/wow_search_20261122/symmetry/README.md#reproduction)

## Selected verification entry points

| Study | What to reproduce |
| --- | --- |
| Mathematical art | [Two orbital exposures and their time ordering](../results/art_20261122/rosette/README.md#reproduce) · [All four studies and their code](../results/art_20261122/README.md) |
| Spatial search of 2026-11-21 | [Regenerate one trajectory or all 30 portrait inputs](shape-search-20261121.md#reproduce-the-portraits) |
| Optical cavity, Newton basins and geodesic flow | [Optical field and basin checks](../results/wow_search_20261122/orbit/README.md#selection-and-replication) · [Fifteen Newton outcomes](../results/wow_search_20261122/boundary/README.md#reproduction) · [Straight periodic ray](../results/wow_search_20261122/lattice/README.md#reproduction) |
| Product-posterior maps | [Fate and singular preimages](../results/product-attractor-status/README.md#storage-and-reproduction) · [Unstable branches and returns](../results/product-larger-attractor/README.md#reproduction) · [Instability and folds](../results/product-mechanism/README.md#reproduction) |
| Icosahedral flow | [Search and absorbing bounds](../research/icosahedral_search/README.md) · [Equilibria](icosahedral-structure.md) · [Numerical returns](../research/icosahedral_returns/README.md) |
| RAAR symbolic dynamics | [Certificate and proof](raar-symbolic-dynamics/README.md) |
| Matrix scaling | [Overrelaxed Sinkhorn](../results/sinkhorn-screen/README.md) · [Unequal-rate inputs and continuations](../results/sinkhorn-unequal-rates/README.md) |
| Adaptive solvers | [RK45](rk45-controller/README.md) · [RK23](rk23-controller/README.md) |
| Acceleration and optimizer memory | [Scalar acceleration](../research/unexpected_search/README.md) · [Yogi](yogi-memory/README.md) · [Adamax](adamax-geometry.md) |
| Motor periodic attractor | [Five-cluster certificate](motor-five-cycle-certificate.md) · [Proof implementation](motor-proof-review.md) |
| Structured-data FastICA | [Algebra and numerical controls](../results/fastica-tripling/README.md) · [Local noninjectivity](../results/fastica-local-injectivity/README.md) · [Rejected second preimage](../results/fastica-preimage-rejection/README.md) |
| Hydrogel filaments | [Planar solver](../results/hydrogel-planar-controls/README.md) · [Spatial solver](../results/hydrogel-spatial-controls/README.md) · [Thin rods](../results/hydrogel-slender-screen/README.md) · [Refinement and tangent checks](../results/hydrogel-variational-checks/README.md) |

The [results index](../results/README.md) lists additional studies. To restore a large input, select its study or filename in the [cache instructions](data-storage.md). Older verifiers may need `--materialize` to restore original paths or [archived figure exports](manuscript-editing.md).

## Recorded environments and runtime checks

The original GPU and CPU controls used Python 3.14.3. From the repository root, create a virtual environment with `python -m venv .venv`. Activate it with `source .venv/bin/activate` on Unix or `.\.venv\Scripts\Activate.ps1` in PowerShell, then install the base packages. GPU controls require an OpenCL driver with double-precision support. [Recorded package versions](../requirements-lock.txt)

```sh
python -m pip install -r requirements.txt
python -m pip check
```

For the archived DCR source solver, also run `python -m pip install -r requirements-dcr.txt`. This installs scikit-sundae 1.2.0 with the SUNDIALS bindings. The CPU check ran for 0.01 time units and matched the reference to within 1.16 × 10⁻¹⁰. [DCR commands and scope](dcr-research.md) · [Dependency check at `3cf58d0`](../records/optional-dcr-packaging.json)

| Revision | Verified |
| --- | --- |
| `8a3564c` | Four GPU/CPU runtime controls in a fresh environment. [Installation and checks](../records/packaging-validation.json) |
| `a514509` | The same four controls, optional DCR control and package compatibility in a fresh public clone. [Validation](../records/public-clone-validation.json) · [File-copy audit](../records/publication-copy-integrity.json) |
| `ad83bf2` | Three- and eight-row motor simulations compared on GPU and CPU in isolated Python. [Validation](../records/motor-public-clone-validation.json) |
| `bcd5119` | All 24 five-cluster proof files matched their hashes. Algebra and exact symmetry checks passed. [Proof export](../records/motor-session-publication-20261030.json) |

## Rebuild figures from saved data

For current main-page images, see the [portrait renderer](../tools/README-portraits.md). The following scripts rebuild individual archived figures or observables in place. Run the script for the figure you want from a disposable checkout.

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

# Reproduce finite-sample FastICA dynamics

Use Python 3.12 and the versions in
[environment.json](../../results/fastica/environment.json). Temporary arrays
and trajectories stay outside Git. The portable scripts change only path
resolution and line endings; exact executed originals and a hash mapping are
in [the evidence packet](../../results/fastica/README.md).

From the repository root, initialize a **new** directory with the five exact
published input matrices and initial vectors:

```powershell
$env:CHAOS_FASTICA_OUTPUT='D:/Codex/fastica-reproduction'
$env:NUMBA_CACHE_DIR='D:/Codex/fastica-numba'
$env:MPLCONFIGDIR='D:/Codex/fastica-matplotlib'
$env:NUMBA_NUM_THREADS='6'
$env:OPENBLAS_NUM_THREADS='1'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONUTF8='1'
python research/fastica/paths.py
python research/fastica/long.py
python research/fastica/precision.py --root $env:CHAOS_FASTICA_OUTPUT --digits 70
python research/fastica/precision.py --root $env:CHAOS_FASTICA_OUTPUT --digits 100 --cases 204 351
python research/fastica/library_probe.py
python research/fastica/controls.py
python research/fastica/stopping.py
python research/fastica/basin.py
python research/fastica/perturb.py
python research/fastica/portrait_data.py
python research/fastica/render_histograms.py
```

The library observation hook is version-sensitive and asserts that it has
found the intended first-unit frame. Use scikit-learn 1.9.1. The underlying
iteration is not patched. Intermediate `rates` in the library receipt are
budget-normalized accumulations; use `complete_window_rates` and
`window_counts` when a run stops early.

To repeat the original 576-case search, use a separate new output directory,
run `pilot.py --out <directory>/pilot-01`, set `CHAOS_FASTICA_OUTPUT` to that
directory, and run `follow.py`. The original selection uses four measured
1,000-step windows, a leading-rate threshold of 0.005, the stated recurrence
test and a minimum update norm of $10^{-8}$. The five selected exact inputs
are retained independently of future screening variability.

`portrait_data.py` generates one million consecutive states for each of four
retained cases. Only Cases 204 and 351 are used in the published portrait.
`render_histograms.py` retains every point in the declared projection; it
writes small count arrays and camera receipts. The gallery uses those counts
through the common `tools/render_portraits.py` renderer. Do not replace the
histograms with thinned trajectory files under their existing names.

Floating-point chaotic trajectories can separate across compilers and library
versions. Compare equations, window statistics and outcome classifications;
byte identity is required for the frozen evidence, not for every fresh run.
The pilot's first step check was exploratory; the archived `controls.py` run
is the reproducible 320-state implementation check.

```text
python research/fastica/verify.py
```

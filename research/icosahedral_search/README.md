# Reproduce the icosahedral flow search

The portable scripts preserve the original calculations and change only path
resolution. Original sources are archived in `results/icosahedral-flow/source-original`.
The shared Gaussian field and full variational integrator come from the frozen
tetrahedral search. Their hashes are recorded as inputs.

Use Python 3.12 with NumPy, SciPy, Numba, Matplotlib, Pillow and mpmath; exact
versions are in the packet's `environment.json`. Select a **new** output folder:

```powershell
$env:CHAOS_ICOSA_OUTPUT='D:/Codex/reproduce-icosahedral'
$env:CHAOS_SEARCH_OUTPUT='D:/Codex/reproduce-icosahedral/shared-cache'
$env:TEMP='D:/Codex/temp'
$env:TMP=$env:TEMP
$env:MPLCONFIGDIR='D:/Codex/matplotlib-cache'
$env:NUMBA_CACHE_DIR='D:/Codex/numba-cache'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONUTF8='1'
New-Item -ItemType Directory -Force -Path $env:TEMP | Out-Null
```

From the repository root, run `python research/icosahedral_search/<name>.py`:

1. `icosa_pilot`, `icosa_gallery`.
2. `icosa_spectra`, `icosa_followup_gallery`, `icosa_long`.
3. `icosa_dop853`, `icosa_dop853_177`, `icosa_basins`, `icosa_cycle`.
4. `icosa_render_long`, `icosa_final_visuals`, `icosa_chirality`.
5. `icosa_algebra`, `icosa_absorbing`.
6. `reiter_replay`, `icosa_compare`.

Stages refuse to replace existing run directories. Allow at least 1 GB of scratch
space. Long integrations use four CPU threads. The geometric comparison also
reads the frozen 108-regime catalog, six earlier candidates and two tetrahedral
clouds already present in the repository.

The replay implements Reiter's mathematical formula independently; its four
published rounded iterates are checked before producing a reference cloud.
Figure 2 is not reproduced because its parameters are unavailable in the appendix.

Verify the committed packet without repeating the search:

```powershell
python research/verify_icosahedral.py --scratch D:/Codex/verify-icosahedral
```

This checks provenance, screen accounting, spectra, failures, rendering geometry
and full point counts, then re-evaluates the absorbing-ball interval inequalities.
It does not certify asymptotic chaos or historical first discovery. Floating-point
orbits can separate across platforms; statistical reproduction matters alongside
source and data hashes.

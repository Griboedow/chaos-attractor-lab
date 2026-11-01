# Reproduce the two-strip proposal

Use the versions in `results/icosahedral-returns/environment.json`. Inputs
default to the three older frozen packets; no long trajectory needs regeneration.
Set an explicit fresh output directory and keep caches on D:

```powershell
$env:CHAOS_ICOSA_RETURNS_OUTPUT='D:/Codex/reproduce-icosahedral-returns'
$env:TEMP='D:/Codex/temp'
$env:TMP=$env:TEMP
$env:MPLCONFIGDIR='D:/Codex/matplotlib-cache'
$env:NUMBA_CACHE_DIR='D:/Codex/numba-cache'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONUTF8='1'
New-Item -ItemType Directory -Force -Path $env:TEMP | Out-Null
```

From the repository root, run in this order:

```text
python research/icosahedral_returns/icosa_return_scout.py
python research/icosahedral_returns/icosa_return_analysis.py
python research/icosahedral_returns/icosa_return_strips.py
python research/icosahedral_returns/icosa_return_controls.py
python research/icosahedral_returns/icosa_return_grid.py
python research/icosahedral_returns/icosa_spherical_prior.py
python research/icosahedral_returns/icosa_return_visuals.py
```

Stages refuse existing run directories. `icosa_poincare.py` and
`icosa_fast_return.py` are imported kernels. Optional input overrides are
`CHAOS_ICOSA_INPUT`, `CHAOS_ICOSA_STRUCTURE_INPUT` and
`CHAOS_ICOSA_CYCLES_INPUT`; see `return_paths.py`.

The spherical comparison reads the existing 80,000-state sample with stride 8.
At packaging time this was checked bitwise against stride 128 of the original
1,280,000-state array, matching the executed scratch source exactly.

For verification without rerunning the complete search:

```text
python research/verify_icosahedral_returns.py --scratch D:/Codex/verify-icosahedral-returns
```

This repeats selected independent numerical controls; it is not a validated ODE
solver or a proof of chaos. Record the date of any new published image with
`tools/render_dates.py record`, then synchronize the dated captions.

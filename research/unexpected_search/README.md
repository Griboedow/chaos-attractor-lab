# Reproduce the tetrahedral-flow and scalar-memory searches

The equations and numerical routines are the same as the archived original scripts. The only portability changes resolve the repository and output directories through `paths.py`. Exact originals, including the D-drive paths used during the investigation, are in `results/unexpected-search/source-original`.

Choose a **new** output directory. Several stages intentionally refuse to overwrite an existing run. On Windows, keep outputs and caches on D:

```powershell
$env:CHAOS_SEARCH_OUTPUT='D:/Codex/reproduce-unexpected-search'
$env:TEMP='D:/Codex/temp'
$env:TMP=$env:TEMP
$env:MPLCONFIGDIR='D:/Codex/matplotlib-cache'
$env:NUMBA_CACHE_DIR='D:/Codex/numba-cache'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONUTF8='1'
New-Item -ItemType Directory -Force -Path $env:TEMP | Out-Null
```

Use Python 3.12 with NumPy, SciPy, Numba, Matplotlib, Pillow and mpmath; exact versions are in the packet's `environment.json`. From the repository root, run these scripts with `python research/unexpected_search/<name>.py` in the following dependency order:

1. `anderson_pilot`, `anderson_extended`, `anderson_gallery`, `anderson_long`, `anderson_spectra`, `anderson_precision`, `anderson_basins`.
2. `tetra_vortex`, `tetra_gallery`, `tetra_long`, `tetra_spectra`, `tetra_dop853`, `tetra_basins`.
3. `long_galleries`, `tetra_render`, `symmetry_audit`, `final_visuals`, `compare_unexpected`, `summarize_evidence`.

The geometric comparison reads the repository's already frozen 108-regime catalog and six earlier candidates. Their exact input hashes are recorded in `reference-inputs.json`. The rendering stage temporarily stores two 3.2-million-point arrays; allow at least 1 GB of working disk space. CPU kernels use four threads. Floating-point trajectories separate chaotically across compiler/library changes, so exact byte reproduction is not assumed across platforms; compare statistics and equations as well as provenance.

To validate the committed evidence without rerunning the searches:

```text
python research/verify_unexpected_search.py
```

This verifies hashes, screen accounting, interval root margins, numerical trace identities, arbitrary-precision results, camera rigidity and actual-state histogram counts. It cannot verify historical first discovery.

# Reproduce the structural follow-up

These scripts change only path resolution from the archived original sources.
They use the same environment as the earlier icosahedral-flow packet. All output
requires an explicit location, preferably on D.

For complete reproduction, first regenerate the dense arrays by running the
earlier [icosahedral search](../icosahedral_search/README.md) through `icosa_long`
and `icosa_final_visuals`, including their dependencies. Point `CHAOS_ICOSA_INPUT`
at that output directory. It must contain `icosa-v1/input.npz`, the dense arrays
under `icosa-long`, and the frozen-camera data `icosa-final/portrait-177.npz`.

```powershell
$env:CHAOS_ICOSA_INPUT='D:/Codex/reproduce-icosahedral'
$env:CHAOS_ICOSA_STRUCTURE_OUTPUT='D:/Codex/reproduce-icosahedral-structure'
$env:CHAOS_SEARCH_OUTPUT='D:/Codex/reproduce-icosahedral-structure/shared-cache'
$env:TEMP='D:/Codex/temp'
$env:TMP=$env:TEMP
$env:MPLCONFIGDIR='D:/Codex/matplotlib-cache'
$env:NUMBA_CACHE_DIR='D:/Codex/numba-cache'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONUTF8='1'
New-Item -ItemType Directory -Force -Path $env:TEMP | Out-Null
```

Run from the repository root, in this dependency order:

```text
python research/icosahedral_structure/icosa_equilibria.py
python research/icosahedral_structure/certify_icosa_equilibria.py --precision 60
python research/icosahedral_structure/certify_icosa_equilibria.py --precision 90
python research/icosahedral_structure/icosa_saddle_starts.py
python research/icosahedral_structure/icosa_itinerary.py
python research/icosahedral_structure/icosa_graph_geometry.py
python research/icosahedral_structure/icosa_mechanism_visuals.py
```

Use a new output directory; most stages refuse an existing run folder. The two
certificate precision levels intentionally write separate files in the inventory
folder. Do not point reproduction output at the frozen packet.

To verify the committed evidence and recertify its root boxes without regenerating
the long orbit:

```text
python research/verify_icosahedral_structure.py --scratch D:/Codex/verify-icosahedral-structure
```

The verifier creates a temporary directory within the explicit scratch location,
copies the archived root centers there, and evaluates all 262 interval boxes anew.
It does not prove the inventory exhaustive or certify the chaotic component.

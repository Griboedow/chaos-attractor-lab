# Reproduce the numerical cycle follow-up

Install the versions in `results/icosahedral-cycles/environment.json`; the new
dependency is `python-flint==0.9.0`. Its use here is numerical, not a validated
ODE certificate. Set all temporary and cache locations explicitly on D.

For the complete recurrence and neighborhood audit, regenerate
`icosa-long/177-dense.npy` using the earlier
[search pipeline](../icosahedral_search/README.md). Other stages read its frozen
coefficient packet and the frozen structural packet without modification.

```powershell
$env:CHAOS_ICOSA_INPUT='D:/Codex/reproduce-icosahedral'
$env:CHAOS_ICOSA_STRUCTURE_INPUT='D:/Codex/chaos-attractor-lab/results/icosahedral-structure'
$env:CHAOS_ICOSA_CYCLES_OUTPUT='D:/Codex/reproduce-icosahedral-cycles'
$env:CHAOS_SEARCH_OUTPUT='D:/Codex/reproduce-icosahedral-cycles/shared-cache'
$env:TEMP='D:/Codex/temp'
$env:TMP=$env:TEMP
$env:MPLCONFIGDIR='D:/Codex/matplotlib-cache'
$env:NUMBA_CACHE_DIR='D:/Codex/numba-cache'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONUTF8='1'
New-Item -ItemType Directory -Force -Path $env:TEMP | Out-Null
```

Run from the repository root in this order:

```text
python research/icosahedral_cycles/icosa_recurrences.py
python research/icosahedral_cycles/icosa_relative_cycles.py
python research/icosahedral_cycles/icosa_cycle_refinement.py
python research/icosahedral_cycles/icosa_itinerary_thresholds.py
python research/icosahedral_cycles/icosa_directed_cycles.py
```

Then run all eight high-precision controls:

```powershell
foreach ($idx in @(0,20)) {
    foreach ($step in @('.04','.02')) {
        python research/icosahedral_cycles/icosa_taylor_precision.py --index $idx --full --order 24 --step $step --digits 80
    }
    foreach ($precision in @(80,110)) {
        python research/icosahedral_cycles/icosa_taylor_precision.py --index $idx --full --order 32 --step .04 --digits $precision
    }
}
python research/icosahedral_cycles/icosa_taylor_checks.py
python research/icosahedral_cycles/icosa_cycle_visuals.py
```

Stages refuse existing run folders or individual precision files. Use a fresh
output location. The original executed scripts and portable scripts differ only
in path resolution. Finite optimizer results can vary slightly across runtimes.

To verify the committed packet without regenerating the long chaotic trajectory:

```text
python research/verify_icosahedral_cycles.py --scratch D:/Codex/verify-icosahedral-cycles
```

This checks stored evidence and repeats local numerical controls; it does not
enclose ODE truncation errors or prove asymptotic chaos or novelty.

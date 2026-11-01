# A Whitney fold in the attracting spatial-circle RAAR set

[Manuscript](manuscript.pdf) · [Current report](../../docs/raar-critical-membership.md) · [Full proof](../../results/raar-critical-membership/theorem.md) · [Frozen packet](../../results/raar-critical-membership/README.md)

The result places a smooth fold inside the previously certified attracting
positive-entropy set. It excludes ambient conjugacy on an open neighborhood
of that set to a map locally injective throughout the corresponding set.
Noninjectivity restricted to K, transitivity, cloud equality and historical
uniqueness remain unproved. The frozen trapping paper is unchanged.

## Verify the archived evidence

Run from the repository root with Python, NumPy, Numba, mpmath and Matplotlib
available. The original runtime used Python 3.12.14, NumPy 2.5.4, Numba 0.68.1
and mpmath 1.3.0. Choose unused output folders.

The verifier also checks the older frozen symbolic-dynamics packet and uses
the original reference trajectory. Restore just these optional historical
inputs at their ignored legacy paths first:

```powershell
python tools/data_cache.py fetch --group raar-symbolic-dynamics --materialize
python tools/data_cache.py fetch results/raar-circles/render/case-1051-animation.npz --materialize
```

The nine optional proposal arrays listed in this packet's `omissions.json`
are not required. The complete final figure input is retained as an explicitly
registered essential input. See [storage instructions](../../docs/data-storage.md).

```powershell
python research/verify_raar_critical_membership.py --out outputs/critical-verification
```

This checks every archived file, replays the complete accounting audits and
numerical figure, and compares analytic/automatic derivative enclosures at
fixed test states. It does not regenerate all proof images. The mathematical
crossing-continuum argument is in the manuscript and theorem.md.

## Recompute every new interval enclosure

The following commands use the frozen exact window. The input path is long
but identical for every calculation.

```powershell
$criticalInput = 'results/raar-critical-membership/runs/raar-critical-thin-window01/input.json'
python research/certify_raar_coverings.py --input $criticalInput --edge 0 --out outputs/critical-replay/raar-critical-thin-cover01
python research/check_raar_coverings_ad.py --input $criticalInput --edge 0 --out outputs/critical-replay/raar-critical-thin-cover-ad01
python research/certify_raar_critical_doubleton.py --input $criticalInput --out outputs/critical-replay/raar-critical-thin-faces01
python research/check_raar_critical_crossing_ad.py --input $criticalInput --out outputs/critical-replay/raar-critical-thin-faces-ad01
python research/certify_raar_critical_crossing.py --input $criticalInput --out outputs/critical-replay/raar-critical-thin-faces-natural01
python research/certify_raar_critical_fold_region.py --input $criticalInput --out outputs/critical-replay/raar-critical-thin-fold-region01
python research/check_raar_critical_fold_ad.py --input $criticalInput --out outputs/critical-replay/raar-critical-thin-fold-ad01
New-Item -ItemType Directory -Path outputs/critical-replay/raar-critical-thin-window01
Copy-Item -LiteralPath $criticalInput -Destination outputs/critical-replay/raar-critical-thin-window01/input.json
python research/audit_raar_critical_membership.py --base outputs/critical-replay --fold raar-critical-thin-fold-region01 --out outputs/critical-replay/audit-primary
python research/audit_raar_critical_membership.py --base outputs/critical-replay --fold raar-critical-thin-fold-ad01 --out outputs/critical-replay/audit-secondary
```

The two map/derivative realizations share the interval library, mean-value
enclosure and partition driver. Numerical QR only chooses finite matrices;
the interval inverse and linear mismatch terms establish inclusion. Different
platform choices may produce different valid subdivisions. A failed or
incomplete run is not a certificate. Recompute the older trapping dependency
using [its own commands](../raar-trapping/README.md).

## Figure and publication

```powershell
python research/render_raar_critical_membership.py --proof results/raar-critical-membership/runs/raar-critical-thin-fold-region01 --out outputs/critical-figure
python research/build_raar_critical_manuscript.py --output outputs/critical-paper/manuscript.pdf
python research/audit_raar_critical_manuscript.py
```

The builder consumes the frozen final figure and evidence. Editable prose is
in `content.py`; review.json records visual inspection of every final page.
The figure is a diagnostic and does not introduce another gallery candidate.
The actual trajectory portrait and slow camera rotation remain on the
[repository main page](../../README.md).

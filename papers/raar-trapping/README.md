# An attracting positive-entropy set in spatial-circle RAAR

[Manuscript](manuscript.md) · [Current report](../../docs/raar-trapping.md) · [Frozen evidence](../../results/raar-trapping/README.md)

This companion to the [symbolic-dynamics study](../raar-symbolic-dynamics/README.md) constructs an explicit strict trapping region for the exact R1051 spatial-circle RAAR map. Its maximal invariant set attracts a neighborhood and has topological entropy at least `log(2)/49`. The proof uses the established RAAR algorithm and covering theory.

The attracting set need not be transitive. Equality with the numerical portrait, a physical invariant measure and historical uniqueness remain unresolved. The existing least-period-7 and least-period-49 saddles are inside this attracting set. The seven-cycle's negative determinant gives only the stated local smooth exclusion of fixed-time ODE maps on R3.

## Restore and audit

The original environment used Python 3.12.14, NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1 and mpmath 1.3.0. Rendering additionally uses Matplotlib. Use a fresh output directory on a disk with at least 1 GB free. Run from the repository root:

```powershell
python research/verify_raar_trapping.py --out outputs/trapping-verification
```

This restores the exact packed mask and sorted keys, checks the complete source-range accounting, rechecks exact window/cycle inclusion, audits the prior symbolic certificates and replays the finite figure data. It does **not** recompute all 21 million outward image enclosures.

## Recompute every enclosure

```powershell
python research/restore_raar_trapping.py --out outputs/trapping-replay/input
python research/audit_raar_trapping_batch_fast.py --input outputs/trapping-replay/input --start 0 --stop 5351450 --out outputs/trapping-replay/part0
python research/audit_raar_trapping_batch_fast.py --input outputs/trapping-replay/input --start 5351450 --stop 10702901 --out outputs/trapping-replay/part1
python research/audit_raar_trapping_batch_fast.py --input outputs/trapping-replay/input --start 10702901 --stop 16054351 --out outputs/trapping-replay/part2
python research/audit_raar_trapping_batch_fast.py --input outputs/trapping-replay/input --start 16054351 --stop 21405802 --out outputs/trapping-replay/part3
python research/audit_raar_trapping_closure.py --input outputs/trapping-replay/input --batches outputs/trapping-replay/part0 outputs/trapping-replay/part1 outputs/trapping-replay/part2 outputs/trapping-replay/part3 --out outputs/trapping-replay/accounting
```

The four batch commands are independent and may run concurrently. Each compiles the outward kernels without cache, records its sources, checks source identity on completion and preserves every completed chunk. Some source ranges are much slower near the curved projection geometry. A run that reports added cells or unresolved geometry is not a closure certificate.

The final raster kernel stops looking for further missing targets as soon as subdivision is already necessary. The original raster control flow is retained, together with all completed common-chunk comparisons. These share the same arithmetic and geometric enclosure; they are not independent proof stacks.

## Figures and manuscript

```powershell
python research/render_raar_trapping.py --input outputs/trapping-replay/input --out outputs/trapping-replay/figures
```

Both new figures were drawn on **2026-11-01**, recorded in the [drawing registry](../../provenance/render-dates.json). The earlier [rotation](../../docs/assets/raar-circles-rotation.gif) shows the native geometry.

The manuscript explicitly states the floating-point arithmetic assumptions, the compactness/attraction argument, the transfer of the covering entropy bound, qualitative parameter persistence and unresolved equivalence questions. The additional primary-source review is in the [priority record](../../results/raar-trapping/priority.json). Third-party publication PDFs are not redistributed.

The complete text is edited directly in [manuscript.md](manuscript.md). Run `python tools/check_manuscripts.py` after editing. The original build, review and manifest records describe the historical publication edition; see [editing and archive instructions](../../docs/manuscript-editing.md). Frozen numerical verifiers may require the explicitly restored files described there.

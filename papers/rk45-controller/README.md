# Twin rosettes and two expanding directions in an adaptive Dormand–Prince solver

**[Research manuscript (Markdown)](manuscript.md)** · [Figures](figures/) · [Experiment report](../../docs/rk45-oscillator.md) · [Priority audit](../../docs/rk45-priority.md)

This companion research manuscript examines selected numerical regimes in the existing SciPy RK45 algorithm applied to a stable linear spring. Its portraits show the actual three-dimensional numerical state: normalized position, normalized velocity and proposed step.

![Twin rosettes and split fans, actual accepted-step states](../../docs/assets/rk45-oscillator-geometry.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

The [rotation](../../docs/assets/rk45-oscillator-rotation.gif) moves only a rigid camera around fixed data. The paper also shows initial volumes approaching the clouds, independent expansion estimates and two declared parameter slices.

## Claims and limits

| Result | Status |
| --- | --- |
| D335 twin rosettes | Two positive finite growth estimates across all selected long-run, independent 70/90-digit and actual-library windows |
| D448 split fans | One positive finite estimate across the corresponding independent checks |
| Initial-state controls | Every selected nonzero solution/step start retains leading expansion; nearby parameters also include periodic outcomes |
| Narrow D335 damping slice | 125 of 257 settings retain two positive finite estimates; 63 are numerically periodic |
| Accepted-solution bounds | Sufficient bound proved, with exact rational certificates for three selected settings |
| Near-half-turn irregular oscillation and adaptive norm dependence | Established prior mechanisms, explicitly credited |
| Full-state trapping, asymptotic chaotic attraction, historical uniqueness | Unresolved |

These are numerical candidates in a known algorithm. The paper does not establish a new solver or a new attractor class. Exact algebra, finite measurements and unresolved priority are distinguished throughout. This is a research draft for author and scientific review; no external peer review or author identities are implied.

## Reproduction

The figures use retained histograms, window estimates and actual-state snapshots. Audit the scientific packet and independently recompute the exact polynomial checks with SymPy:

```bash
python research/verify_rk45_controller.py --root results/rk45-controller --assets docs/assets
python research/audit_rk45_structure.py --cases results/rk45-controller/source/structure/input-cases.json --out output/rk45-structure-audit.json
```

The complete text is edited directly in [manuscript.md](manuscript.md). Run `python tools/check_manuscripts.py` after editing. The original build, review and manifest records describe the historical publication edition; see [editing and archive instructions](../../docs/manuscript-editing.md). Frozen numerical verifiers may require the explicitly restored files described there.

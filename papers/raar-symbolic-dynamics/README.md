# Certified symbolic dynamics in a spatial-circle RAAR map

[Manuscript](manuscript.pdf) · [Research report](../../docs/raar-symbolic-dynamics.md) · [Frozen packet](../../results/raar-symbolic-dynamics/README.md)

Seven pages, three figures, editable prose in `content.py`. The result is a compact invariant positive-entropy set in an exact existing RAAR map. Attraction, equality with the numerical cloud and historical uniqueness remain unresolved.

## Results and evidence

- Eight covering relations for G=F^7, independently realized with vector/analytic and scalar/automatic interval engines: 16,467 and 9,449 accepted cells.
- Complete dyadic tilings, strict exit/entry conditions, regularity, invertibility and N0/N1 disjointness.
- Arbitrary binary routes with 49 original iterations per symbol; topological entropy at least log(2)/49 > 0.01414.
- A separate located-root certificate for a least-period-49 saddle, with exactly one unstable and two stable multipliers.
- All proposals, the interrupted pilot and all figure revisions retained. Long numerical monodromy products are not used to classify stability.

The graph-to-entropy argument uses established covering theory, cited in the paper. It does not certify an attracting chaotic set or a new equivalence class.

## Reproduction

Run from the repository root. Python 3.12, NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1 and mpmath 1.3.0 were used for numerical and interval work. The full interval scripts need only mpmath; the output auditor needs NumPy and SciPy. PDF rebuilding additionally uses Matplotlib, ReportLab and pypdf.

```powershell
python research/verify_raar_symbolic.py --output symbolic-audit.json
python research/audit_raar_symbolic_manuscript.py
```

To recompute every interval enclosure, choose a fresh output directory. These commands create subdirectories and refuse to overwrite existing runs:

```powershell
$windowInput = 'results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json'
foreach ($edge in 0..7) {
    python research/certify_raar_coverings.py --input $windowInput --edge $edge --digits 70 --out "outputs/symbolic-replay/primary-$edge"
    python research/check_raar_coverings_ad.py --input $windowInput --edge $edge --digits 100 --out "outputs/symbolic-replay/second-$edge"
}
python research/certify_raar_cycles.py --input results/raar-symbolic-dynamics/runs/raar-symbolic-orbits01/certificate-input.json --out outputs/symbolic-replay/cycle
python research/audit_raar_symbolic_cycle.py --windows $windowInput --cycle outputs/symbolic-replay/cycle/symbolic-period-49.json --out outputs/symbolic-replay/cycle-audit
```

All engine certificates preserve their local source snapshots and input hashes. The second engine shares the subdivision driver and interval library with the first; they are not fully independent proof stacks.

To rebuild the manuscript and figures:

```powershell
python research/render_raar_symbolic.py --base results/raar-symbolic-dynamics/runs --out outputs/symbolic-replay/figures
python research/build_raar_symbolic_manuscript.py --output outputs/symbolic-replay/paper/manuscript.pdf
```

The manuscript builder consumes the frozen final figures. The renderer reconstructs all one million density contributions, exact plotting coordinates and distances from retained arrays. Rendering changes do not modify the evidence packet. `build.json`, `review.json` and `manifest.json` identify the reviewed publication bytes.

Priority review: [main audit](../../docs/raar-priority.md) and [new method/source record](../../results/raar-symbolic-dynamics/priority.json). No exact matching covering result was located in the bounded search; that is not a priority proof.

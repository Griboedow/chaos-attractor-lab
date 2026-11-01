# Tolerance-scale geometry and numerical hyperchaos in an adaptive Runge-Kutta solver

**[Research manuscript (PDF)](manuscript.pdf)** · [Editable text](content.py) · [Figures](figures/) · [Detailed experiment report](../../docs/rk23-oscillator.md) · [Priority audit](../../docs/rk23-priority.md)

The 11-page research manuscript studies intricate accepted-step dynamics in the **existing standard RK23 solver** applied to a stable linear spring. The physical solution decays. The displayed geometry belongs to the numerical algorithm's position, velocity and proposed step.

![Actual RK23 accepted-step states](figures/figure-1-geometry.png)

The [rotating view](../../docs/assets/rk23-oscillator-rotation.gif) changes only the camera. All plotted points are actual recorded states, with equal spatial units. The paper also shows how 32,768 initial states per case approach the sampled clouds.

## Claims and limits

| Result | Status |
| --- | --- |
| Explicit bound on accepted solution values | Proved under a sufficient relative-tolerance condition, independent of the proposal policy |
| Sharper bound for O526 | Exact-rational Sturm and gain certificate, with the induction written out |
| Folded-triangle and scissor-wing regimes | Five-million-step evidence, independent 70/90-digit checks and actual SciPy execution |
| Two expanding directions in O526 | Positive finite estimates in every declared independent window |
| O37 and O689 second directions | Positive long-window estimates, but mixed shorter independent windows; not promoted under the stricter criterion |
| Full-state trapping region, asymptotic chaotic attraction or historical uniqueness | Unresolved |

Adaptive-controller instability, tolerance-scale oscillation and error-norm dependence are established prior work. The manuscript credits that literature and records full-text access gaps. Its proposed contribution is a concrete, reproducible description of selected regimes and accompanying solution bounds; it does not establish a new attractor class. No author identities or affiliations were inferred.

## Reproduction

The build reads the frozen results/rk23-controller packet; it does not rerun a trajectory to obtain the figures. Plot lines and labels are vector PDF objects. Density portraits and dense approach clouds are raster layers derived from the archived data.

    python research/build_rk23_manuscript.py --output output/pdf/rk23-controller-manuscript.pdf
    python research/audit_rk23_manuscript.py
    python research/verify_rk23_controller.py --root results/rk23-controller --assets docs/assets

The first command creates a new review artifact and figure directory beside it. The latter commands audit the publication copy and frozen scientific packet. Building requires the numerical environment plus ReportLab and pypdf. The independent exact-symbolic check uses SymPy:

    python research/audit_rk23_structure.py --out output/rk23-structure-audit.json

Every manuscript page was rendered and visually inspected. The [review record](review.json) identifies the checked PDF and rendered-page hashes. Numerical claims were compared with the retained arrays, including weaker outcomes. The manifest covers the PDF, editable text, figures and build/review metadata. This is a research draft for author and scientific review, not an externally peer-reviewed publication.

# Yogi: prior work, equivalent descriptions, and defensible claims

Audit date: **1 November 2026**. This is a bounded literature and equation audit, not a certificate of worldwide priority. The [main report](yogi-geometry.md) presents numerical evidence; the [structure note](yogi-structure.md) separates exact identities from observations.

## Claim boundaries

| Possible claim | Status |
|---|---|
| A new optimizer or new update equation | False: Yogi was published in 2018. |
| First chaos in optimization of a scalar quadratic | False: Adam/RMSProp literature already covers this setting. |
| New because the portrait resembles fans or a sunburst | Unsupported: outlines depend on coordinates and camera. |
| First hyperchaos in a three-dimensional discontinuous map | Not claimed; this is an established research area. |
| Persistent two-positive-direction evidence in this specific Yogi recurrence | Supported numerically here, with the declared finite tests and negative outcomes. |
| Previously undescribed Yogi regime | A candidate contribution: no matching report was identified in the inspected material, but coverage is incomplete. |
| A new topological class, inequivalent to every known attractor | Unresolved; no such proof is provided. |

The strongest current article angle is the particular Yogi memory mechanism, its two-direction numerical regime, and a reproducibility failure at a discontinuous initialization. A rounded parameter tuple, fresh name, or new application label is not sufficient novelty.

The [subsequent grid study](yogi-parameter-plane.md) broadens the finite evidence and adds a feathered-arc portrait. A separate [scalar boundedness proof](scalar-optimizer-boundedness.md) covers all the stated memory parameters and nonzero initial memories; its prior-art comparison does not claim that the bounding technique is historically new.

## Original algorithm and convergence results

[Zaheer, Reddi, Sachan, Kale and Kumar (2018), *Adaptive Methods for Nonconvex Optimization*](https://proceedings.neurips.cc/paper/2018/file/90365351ccc7437a1309dc64e4db32a3-Paper.pdf): Algorithm 2 and Section 3.2 were inspected; the [author's 18-page version](https://www.satyenkale.com/papers/yogi.pdf) additionally supplies Theorem 2 and Appendix B. The second-moment lower bound `v_t>=beta2*v_(t-1)` and a bounded-gradient upper bound are already there. Their displayed analysis uses zero first-moment coefficient, globally bounded gradients, and restrictive step/epsilon assumptions. Our fixed-step examples do not satisfy those hypotheses and are not counterexamples to that theorem.

[Ibragimov and Jentzen (2026), *Unified convergence analysis for gradient descent optimization methods in the training of deep neural networks*, arXiv 2607.04233v1](https://arxiv.org/abs/2607.04233v1): Corollary 4.12 explicitly covers Yogi. The full statement was read: its learning rate decays comparably to `n^-nu`, with `3/4<nu<=1`, alongside boundedness and regularity assumptions. A constant learning rate lies outside that statement. “Yogi converges” without its hypotheses would misrepresent this work.

[Dereich, Graeber, Jentzen and Riekert (2025), *Asymptotic stability properties and a priori bounds for Adam and other gradient descent optimization methods*, arXiv 2509.10476v1](https://arxiv.org/abs/2509.10476v1): Theorem 1.3 and Proposition 2.8 were inspected. The related boundedness technique uses an exponential second moment or a dominating preconditioner and a condition corresponding to `beta2>beta1^2`. Yogi's signed accumulator does not generally dominate that exponential average, so direct application requires an additional argument. Our elementary scalar bound is not claimed to originate the general technique.

[Fong and Yang (2026), *Adam at the Edge of Stability: Adaptive Feedback, Provable Oscillation, and Gradient Reversal*, arXiv 2608.20638v2](https://arxiv.org/abs/2608.20638v2) was checked at its revised version, including Lemma H.3 and its proof: its diagonal-quadratic bound requires `beta1^2<beta2`. [Dereich, Do and Jentzen (2026), arXiv 2603.18899v2](https://arxiv.org/abs/2603.18899v2), Corollary 2.7, requires that same strict memory ordering and nonincreasing learning rates tending to zero. Neither statement directly covers the constant-step high-first-memory Yogi examples. [Dimitrieski et al. (2026), arXiv 2603.15823v2](https://arxiv.org/abs/2603.15823v2) already establishes bounded-step robustness for no-momentum RMSProp under strong convexity and global smoothness. These primary results are credited in the boundedness note; [inspection metadata](../results/yogi-plane/bound-sources/) records the exact downloaded versions.

[Dereich, Graeber and Jentzen (2024), *Non-convergence of Adam and other adaptive stochastic gradient descent optimization methods for non-vanishing learning rates*, arXiv 2407.08100v1](https://arxiv.org/abs/2407.08100v1): Theorem 1.2 was read. Its quadratic nonconvergence example has nondegenerate stochastic variability. That result is relevant context, but it does not by itself establish deterministic chaos in the noise-free scalar system studied here.

## Optimizer dynamics already reported

[Bock and Weiss (2022), arXiv 2210.02070](https://arxiv.org/html/2210.02070): Section 4 explicitly reports chaotic behavior of Adam on a scalar quadratic, including fountain-like and Lorenz-like portraits. That directly defeats a broad claim that our unexpected setting alone is new.

[Loh, Mireles James, Xu and Tao (2026), arXiv 2610.30078v1](https://arxiv.org/abs/2610.30078v1), submitted 21 October2026: the introduction, related work, trapping discussion, conclusion, and Adam extension were inspected during the preceding [Adamax audit](adamax-priority.md). It studies nonlinear optimizer dynamics on simple objectives and is essential recent context. A full-text search did not locate Yogi, Adamax, or hyperchaos terminology, but a missing keyword does not exclude equivalent results elsewhere in the paper or literature.

These papers are credited as predecessors. The present study does not rebrand their general phenomenon as a discovery.

[Sebastian Bock's dissertation (repository year 2024), *3D Pose Estimation Problem and a Convergence Analysis of Adaptive Gradient Descent Algorithms*](https://epub.uni-regensburg.de/58326/) was checked separately, including selected sections of Chapters 4–5, the AdaDelta cycle analysis, and the full-state Adam figures 31 and 34. It explicitly omits AdaMax from Chapter 5's analysis. No Yogi or hyperchaos term was located in the extracted text. The broader cycles and chaotic portraits are prior work; the keyword result is not a non-equivalence proof. [Inspection record](../results/optimizer-priority-followup/bock-2024.json).

## Equivalent-language checks

The searches covered not only `Yogi + chaos / attractor / Lyapunov / hyperchaos / bifurcation / limit cycle / quadratic`, but also signed adaptive accumulators, discontinuous maps, border collisions, lines of fixed points, memory maps, and embeddings of classical maps. Some broad Yogi searches mainly returned authors whose surname is Yogi; those are not optimizer evidence.

[Di Marco, Forti, Pancioni, Innocenti and Tesi (2024), *Embedding classic chaotic maps in simple discrete-time memristor circuits*, arXiv 2408.16352v1](https://arxiv.org/abs/2408.16352v1): the primary preprint's introduction and Section 4, including equations 66–73 and Remarks 6–8, were inspected. It shows how familiar logistic/tent and Henon/Lozi dynamics can occur inside higher-dimensional memory circuits. Its three-dimensional example 69 has an explicit first integral and embeds Henon dynamics on invariant level sets. This is a concrete warning against equating a memory variable, a fixed line, or an unexpected application with dynamical novelty.

The full Yogi map is discontinuous and has overlapping inverse branches, whereas those particular Henon/Lozi embeddings are continuous. This excludes a global homeomorphic conjugacy on the full domains. It does not exclude a conjugacy between selected invariant subsets, a factor map, or a different nonsmooth memory construction. No such invariant-set comparison has been completed.

[Linero Bas and Nieves Roldan (2022), arXiv 2209.03582](https://arxiv.org/abs/2209.03582) was inspected in the Adamax audit, including Sections 2/5 and the relevant figures. Its logarithmic and max-type equivalences reinforce the need to compare recurrences after coordinate changes. Merely seeing a `sign`, `max`, or square root does not establish a distinct class.

[Patra, Gupta and Banerjee (2021), *Local and global bifurcations in 3D piecewise smooth discontinuous maps*, DOI 10.1063/5.0010887](https://doi.org/10.1063/5.0010887) is an outstanding close comparator. The [authors' institutional publication listing](https://web.iitm.ac.in/ccsd/page8.html) confirms the paper in this directly relevant class. The publisher returned 403 and the alternative indexed abstract page presented a browser challenge. **Its full equations were not read, and no non-equivalence claim is made against it.** This access gap remains material to a priority claim.

## Implementation identity is also part of the audit

[Optax's Yogi documentation](https://optax.readthedocs.io/en/latest/api/generated/optax.yogi.html) and [source](https://github.com/google-deepmind/optax/blob/main/optax/_src/transform.py) were inspected. The installed 0.2.8 function, its initial state, version metadata, and source hash are recorded separately from the inspected main branch. They differ in first-moment initialization. The actual-library experiment retains its continuous bias correction, rather than switching implementations after warmup. That distinction changes basin selection for case 44 and therefore belongs in the scientific result.

All primary PDFs read locally are recorded by URL and checksum. They are not redistributed. Source snapshots, exact model conventions, every declared numerical outcome, and the initial failed implementation check are retained in the [packet](../results/yogi-geometry/README.txt).

## Wording suitable for the present evidence

“We report reproducible numerical regimes of scalar-quadratic Yogi with intricate full-state geometry. One regime exhibits two positive finite-time branch growth rates across five million iterations, corroborated by independent high-precision finite-separation calculations and an actual Optax run. We also identify an initialization ambiguity at the map's discontinuity. The algorithm and the broader phenomenon of optimizer chaos are established; invariant-set classification, rigorous chaotic attraction, and historical priority for this particular regime remain open.”

This wording can be strengthened only after the missing mathematical and bibliographic work, not by increasing visual polish or renaming the shapes.

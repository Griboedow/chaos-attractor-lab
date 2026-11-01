# Bounded memory and numerical hyperchaos in scalar quadratic optimization

**[Research manuscript (Markdown)](manuscript.md)** · [Figures](figures) · [Detailed experiment report](../yogi-tanh.md) · [Priority audit](../yogi-priority.md)

The manuscript combines a complete scalar boundedness proof with finite numerical evidence for two expanding directions in the **existing** smooth Yogi optimizer. It derives an exact local cycle-stability interval and explicitly identifies the known Adam cycle's local Jacobian equivalence.

![Actual smooth-Yogi states in two rigid camera views](figures/figure-1-geometry.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

The left ribbon is initialization-dependent: its standard library cold start settles. The right crescent is reproduced by actual TensorFlow. The [rotation](../assets/yogi-tanh-rotation.gif) changes only the camera.

## What is proved, and what remains numerical

| Result | Status |
| --- | --- |
| Scalar boundedness under a positive denominator floor and d >= sigma abs(x) | Complete first-exit proof, with no ordering of momentum and second-memory parameters |
| Every finite legal smooth-Yogi orbit is bounded, including bias correction | Consequence of that theorem with an explicit denominator constant |
| Critical fold sheets and global noninjectivity | Exact algebraic consequences of the known map |
| Strict exponential stability interval of the symmetric two-cycle | Proved; its equivalence to the known Adam local Jacobian is explicit |
| Two expanding directions in S0039 and S7239 | Five-million-step numerical checks, independent 70/90-digit calculations, actual unmodified library execution and fresh-start cohorts |
| Asymptotic chaotic attraction or a historically new attractor class | Unresolved |

The algorithm, optimizer chaos on scalar quadratics, and the symmetric-cycle family are prior work. The manuscript does not claim priority of the bounding technique. Its references, inspected source versions and remaining coverage gaps are explicit. No author identities or affiliations were inferred.

## Reproduction

The figures use frozen histograms and numerical summaries from `results/yogi-tanh`. The scientific verifier is:

```bash
python research/verify_yogi_tanh.py --root results/yogi-tanh
```

The complete text is edited directly in [manuscript.md](manuscript.md). Run `python tools/check_manuscripts.py` after editing. The original build, review and manifest records describe the historical publication edition; see [editing and archive instructions](../manuscript-editing.md). Frozen numerical verifiers may require the explicitly restored files described there.

This is a research draft for scientific and author review, not an accepted or externally peer-reviewed publication.

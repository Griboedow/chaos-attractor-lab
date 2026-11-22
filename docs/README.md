# Documentation

The [main page](../README.md) pairs each portrait with its equations and parameters. The [research index](../README.md#research-index) also covers searches that ended in periodic motion, numerical artifacts or open questions.

| Start with | What you will find |
| --- | --- |
| [Four studies from 2026-11-22](chaotic-forms-20261122.md) | A point alternating between rings, a fan in an optical cavity, fifteen Newton basins with one boundary, and a straight orbit on a chaotic surface. |
| [Spatial search from 2026-11-21](shape-search-20261121.md) | Knotted periodic orbits, linked attractors and maps with several connected components. |
| [Reproduce a result](reproduction.md) | A first certificate to check, instructions for individual studies and the software used. |

## Articles, proofs and data

RAAR stands for relaxed averaged alternating reflections. The first six rows cover different properties of the same spatial-circle projection algorithm.

| Study | Report | Full article | Evidence and reproduction |
| --- | --- | --- | --- |
| Spatial-circle RAAR geometry | [Report](raar-circles.md) | [Article](raar-circles/manuscript.md) | [Study notes](raar-circles/README.md) |
| RAAR periodic structure | [Report](raar-cycle-structure.md) | [Article](raar-cycle-structure/manuscript.md) | [Study notes](raar-cycle-structure/README.md) |
| RAAR folding and preimages | [Report](raar-folding.md) | [Article](raar-folding/manuscript.md) | [Study notes](raar-folding/README.md) |
| RAAR symbolic dynamics | [Report](raar-symbolic-dynamics.md) | [Article](raar-symbolic-dynamics/manuscript.md) | [Study notes](raar-symbolic-dynamics/README.md) |
| RAAR attracting set | [Report](raar-trapping.md) | [Article](raar-trapping/manuscript.md) | [Study notes](raar-trapping/README.md) |
| Fold in the RAAR attracting set | [Report](raar-critical-membership.md) | [Article](raar-critical-membership/manuscript.md) | [Study notes](raar-critical-membership/README.md) |
| Motor-cycle replication instability | [Report](motor-replication-index.md) | [Article](replication-instability/manuscript.md) | [Study notes](replication-instability/README.md) |
| RK23 controller dynamics | [Report](rk23-oscillator.md) | [Article](rk23-controller/manuscript.md) | [Study notes](rk23-controller/README.md) |
| RK45 controller dynamics | [Report](rk45-oscillator.md) | [Article](rk45-controller/manuscript.md) | [Study notes](rk45-controller/README.md) |
| Scalar optimizer memory | [Report](yogi-tanh.md) | [Article](yogi-memory/manuscript.md) | [Study notes](yogi-memory/README.md) |
| Gaussian-posterior reflected dynamics | [Report](soft-symbolic-dynamics.md) | [Article](soft-symbolic-dynamics/manuscript.md) | [Study notes](soft-symbolic-dynamics/README.md) |
| Gaussian posterior mean and a sphere | [Main portraits](../README.md#gaussian-posterior-mean-and-a-sphere) | [Numerical article](soft-sphere/manuscript.md) | [Evidence and reproduction](../results/soft-sphere/README.md) |
| Canonical posterior maps | [Main portraits](../README.md#canonical-posterior-maps-ribbon-and-spatial-folds) | [Article and proofs](canonical-sphere/manuscript.md) | [Numerical campaign](../results/canonical-sphere/README.md), [collision and homoclinic certificates](../results/canonical-collision/README.md) |
| Product-alphabet posterior maps | [Main portraits](../README.md#product-alphabets-folded-wings-and-layered-arches) | [Article, entropy, orientation and instability theorems](product-alphabet/manuscript.md) | [Numerical campaign](../results/product-alphabet/README.md), [entropy](../results/product-symbolic/README.md), [orientation and equivalence](../results/product-orientation/README.md), [long horizon and excursions](../results/product-horizon/README.md), [instability and mechanism](../results/product-mechanism/README.md) |
| Nonuniform product posteriors | [Folded-shell portraits](../README.md#folded-shells-from-nonuniform-quantizer-levels) | [Model, controls and prior art](product-alphabet/manuscript.md#32-nonuniform-alphabets-and-folded-shells) | [Complete higher-order and nonuniform searches](../results/product-highorder/README.md) |

## Studies with shorter reports

- **Structured-data FastICA.** A constructed data set gives an attracting set with angle tripling as a factor. The follow-ups distinguish the map near this set from a solenoid embedding and show why a proposed second inverse branch falls outside the set. [Proof and numerical tests](fastica-tripling.md) · [Geometry](fastica-tripling-structure.md) · [Rejected branch](fastica-preimage-rejection.md)
- **Hydrogel filaments.** An independent solver reproduced the planar controls. Spatial runs exposed a coarse-mesh artifact and two close-contact outcomes in thinner rods. Refinement revealed slow drift behind apparent growth. [Eight planar cases](hydrogel-planar-controls.md) · [Spatial controls](hydrogel-spatial-controls.md) · [Thinner rods](hydrogel-spatial-controls.md#circular-section-symmetry-and-thinner-filaments) · [Drift and refinement](hydrogel-variational-checks.md)
- **Sinkhorn matrix scaling.** Long continuations of a 576-input screen exposed transient growth. Of 4,864 rectangular and unequal-rate inputs, two candidates retained irregular motion over 106 million steps. Chaotic attraction remains open. [First screen and absorbing bounds](sinkhorn-screen.md) · [Unequal-rate results and precision tests](sinkhorn-unequal-rates.md)

[Numerical data and cache](data-storage.md) · [Sources and checks](../RECORDS.md) · [Editing articles](manuscript-editing.md)

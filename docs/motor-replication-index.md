# Exact stability index of the replicated five-cluster cycle

The [research manuscript](../papers/replication-instability/README.md) strengthens the earlier [replication limit](motor-five-replication-limit.md). For every integer `k >= 1`, the cycle embedded in `N=8k` rows has exactly

```math
(u,n,s)=(2k-2,1,22k).
```

All stable multipliers satisfy `|mu|<0.13`. Each formerly singleton contrast block has precisely one real multiplier in `(-23.084,-23.083)` and two multipliers within radius `0.13`. Thus the original `N=8` orbit attracts locally, while every enlarged replica is a periodic saddle. The parameters and cluster weights remain fixed. The inherited positive period is certified; its least-period status is unchanged.

![Replication mechanism and the exact unstable count](../papers/replication-instability/figures/figure-3-replication.png)

<!-- drawing-date -->
*First recorded drawing: **2026-10-30** ([record](https://github.com/Griboedow/chaos-attractor-lab/blob/b5afd472f26e7fa247ab0b5dd8644d9ffaabacc6/papers/replication-instability/figures/figure-3-replication.png)).*
<!-- /drawing-date -->

## New proof step

The old sign bracket established at least one unstable root in each new block. The new [spectral checker](../research/certify_motor_replication_index.py) narrows that real-root bracket and encloses synthetic division of the interval characteristic polynomial by the corresponding exact linear factor. Three strict real quadratic Jury inequalities place **both remaining roots** within radius `0.13`. This turns a lower bound into a complete block index.

The calculation is repeated from two separately integrated matrix enclosures: the preserved 512-bit record and a new [640-bit, 192-step, order-54 integration](../results/motor-replication-article/reintegration-640.json). Both [768-bit replays](../results/motor-replication-article/exact-index-recheck.json) give the same conclusion. No interval radius is dropped, and floating-point eigenvalues do not enter the proof.

The tangent space is the direct sum of the 14-dimensional quotient and `km_a-1` copies of the three-dimensional internal block for each group. Counting its invariant blocks yields `13+9(2k-1)+4(k-1)=22k` stable directions, `2(k-1)` unstable directions and the single flow direction. The strict spectral separation and invertible shooting Jacobian also imply persistence of the index in a common sufficiently small parameter neighborhood; no numerical radius of that neighborhood is supplied.

## Nonlinear interpretation

All [eight predeclared N=16 starts](../results/motor-replication-article/nonlinear/protocol.json) are retained: both new blocks, both signs, and perturbation norms `1e-6` and `1e-8`. Full47D integrations reproduce amplification and sign reversal. At the first three returns, all smaller-amplitude predictions agree to relative error below `1e-6`. Tighter DOP853 histories and an independent Radau control pass. These numerical illustrations do not classify the eventual destinations or establish a nonlinear basin bound.

Earlier proof files, source snapshots and raw data are unchanged. The manuscript's [claim map and reproduction guide](../papers/replication-instability/README.md), [literature positioning](../papers/replication-instability/priority-and-positioning.md), vector figures and rotating GIF accompany the new evidence.

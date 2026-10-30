# A certified limit to replication of the five-cluster attractor

The [attracting eight-row cycle](motor-five-cycle-certificate.md) does **not** remain an attracting periodic orbit under exact replication to `N=8k` rows when `k>=2`, at the same normalized parameters and cluster weights. Each of its two singleton cluster states supplies a newly available transverse block with a rigorously enclosed real multiplier in `(-32,-1)`.

This does not change the fixed N=8 theorem. A singleton has no internal contrast direction when it contains just one row. Replication creates those additional directions.

The claim concerns the mathematical family with fixed `rho=129/55`, `r=1+(19/20)/sqrt(32001/500000)`, `b=1` and weights `(2,2,1,2,1)/8`. It does not assert that physical constitutive parameters remain unchanged when the number of motor rows changes.

## Exact conclusion

Replicate the five cluster states with multiplicities

```math
(2k,2k,k,2k,k),\qquad N=8k.
```

Their normalized weights are unchanged, so the certified fourteen-dimensional quotient orbit remains an exact solution. Every cluster of size m has `m-1` independent internal contrast copies. The three original paired-cluster blocks retain their certified stability. Each formerly singleton cluster now has `k-1` copies of its previously absent three-dimensional block.

Both new blocks have a real multiplier strictly between -32 and -1. Consequently the replicated orbit has **at least `2(k-1)` expanding directions** for every integer `k>=2`. This is a lower bound, not a complete spectral count for the new blocks. In particular, the N=16 replica is transversely unstable. No basin, chaotic-attractor or nonlinear bifurcation conclusion follows.

## Validated calculation

The separate [extension result](../results/motor-five-replication-proof/singleton-blocks.json) reads the immutable successful recheck of the N=8 cycle and retains its input hash. It does not repeat or alter the periodic shooting proof. A new subclass initializes the existing weighted field without its old variational blocks, then appends only the two singleton contrast systems. No module globals or frozen proof files are changed.

Three segment integrations use 512-bit Arb arithmetic, 160 subdivisions per segment and Taylor order 50. The original validated shooting boxes contain the exact segment states and duration. If B2 and B4 are the resulting matrices over one relative interval L, the certified permutation swaps these two clusters. Their full returns over the inherited positive period `T=6L` are therefore

```math
M_2=(B_4B_2)^3,\qquad M_4=(B_2B_4)^3.
```

All matrix uncertainties remain in these ordered products. The strict unit-disk Schur tests fail, which alone would be inconclusive. A **separate** characteristic-polynomial sign test gives, for each block,

```math
p(-32)<0,\qquad p(-1)>0.
```

The displayed values are approximately `p(-32)=-9130.5264376294` and `p(-1)=22.08338187555`, with interval radii far smaller than either sign margin. The intermediate value theorem supplies a real root in `(-32,-1)`. This strict root bracket, rather than the failed Schur test, proves transverse instability.

Before integration, 1,692 independently assembled components compare both contrast embeddings with the full 47-dimensional N=16 Jacobian at rational test states. The controls include the cluster containing the eliminated last X coordinate; every component passes. Read-only review also checked the weighted embeddings, new block entries and full return products. The run took approximately 65 seconds.

A separate [768-bit saved-matrix replay](../results/motor-five-replication-proof/matrix-replay.json) independently recomputes both characteristic polynomials and endpoint signs from the archived full matrix enclosures. Both strict root brackets pass again. This is a serialization and spectral replay, not a second integration. Its [checker](../research/check_motor_five_replication_record.py) and executed source are retained alongside the original result.

Reproduce with an unused output path:

```bash
python research/check_motor_five_replication.py --input results/motor-five-cycle-proof/recheck-01.json --output outputs/five-singleton-blocks.json
```

The [implementation](../research/check_motor_five_replication.py), executed source snapshots and [packet manifest](../results/motor-five-replication-proof/manifest.json) preserve the calculation and its provenance. This extension uses the same Picard/Taylor enclosure method and Python-FLINT/Arb arithmetic as its input certificate. It makes no least-period assertion and leaves every earlier packet unchanged.

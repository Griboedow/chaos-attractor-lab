# Approaches to the nonuniform posterior shells

This companion to the [alphabet-extension packet](../product-highorder/README.md)
shows the first 2,000 updates from the same 33 stored initial states for
nonuniform Cases 6 and 7. It does not modify the frozen extension packet.
See the [article and figure](../../docs/product-alphabet/manuscript.md#32-nonuniform-alphabets-and-folded-shells).

The reference consists of every twentieth state from each two-million-state
portrait, giving 100,000 points. Distances are nearest-neighbor Euclidean
distances divided by the reference bounding-box diagonal. Stored arrays
retain all three distance quantiles at every update and the exact plotted
states at the declared logarithmically spaced indices. They do not replace
the full trajectory under its original name. Intermediate full approaches
were held in memory and can be regenerated from the retained inputs.

Independent squared-distance posterior weights and nested reflections
agree with the production formula at both selected steps for all 66
trajectories: 132 checks, maximum absolute error 3.34e-15. This checks the
expression at those states, not accumulated high-precision orbit agreement.
At update 2,000, median relative distances are 0.0006974 and 0.001557.
Finite-cloud distance does not establish an invariant set, basin volume,
or asymptotic chaotic attraction. The plot was visually inspected.

`report.json` records input and output hashes. `manifest.json` covers all
packet files except itself. The existing source-state files can be restored
or regenerated using the alphabet-extension instructions. Reproduce with:

```powershell
python research/render_irregular_approach.py --source D:/runs/irregular-dense --cohort results/product-highorder/runs/product-irregular-cohort01/protocol.json --out D:/runs/irregular-approach --figure docs/assets/diagnostics/product-irregular-approach.png
```

Use a new output directory. The retained Python snapshots fix the executed
recurrence and renderer; the main source tree exposes the same code.

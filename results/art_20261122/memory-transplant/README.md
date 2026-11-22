# Transplanting memory into an optical cavity

Take the instantaneous light field from one run of the cavity model and its slow memory from another. In all 256 pairings of 16 snapshots, the observed detuning regime followed the memory donor. At intermediate starting memories, field U04 changed regime seven times across 21 sampled values. Replays at 640, 1,280 and 2,560 bits gave the same sequence.

![Sixteen optical fields scanned across initial memory. Red cells mark where the implementations disagree on the outcome.](memory-transplant.png)

Follow a row from left to right. Each cell is a fresh run with a larger starting memory. Teal means the memory subsequently stayed low in both observation windows, amber means it stayed high. The right-hand panel enlarges the boundary. The colors and counts are from two binary64 calculations. Red marks their disagreements. At some unmarked cells, both calculations give the same answer and higher precision changes it.

## Experiment

The driven optical cavity map is

```math
u=\frac{|E|^2}{1+|E|^2},\qquad E^+=P+aE\exp[i(\delta+ku+T)],\qquad T^+=(1-h)T+hgu.
```

The squared magnitude of the complex field E is the light intensity. T is the slowly relaxing detuning. In Case 86780 the parameters (P, a, δ, k, g, h) are approximately (0.68380272, 0.47694415, 0.64927274, −17.94019890, 19.38252068, 0.00387439). The model follows the optical ring-cavity feedback studied by Ikeda. [Exact parameters](../../wow_search_20261122/orbit/case-86780.json) · [Model and earlier basin scan](../../wow_search_20261122/orbit/README.md) · [Ikeda (1979)](https://doi.org/10.1016/0030-4018(79)90090-7)

1. Start one trajectory in each observed regime. After 32,000 updates, take eight snapshots from each at 313-step intervals. L01–L08 are the lower-regime donors and U01–U08 are the upper-regime donors. [Starts, step numbers and snapshot coordinates](donors.json)
2. Form all 16 × 16 combinations of a snapshot's field and a snapshot's memory. For a second experiment keep each of the 16 fields and scan 193 evenly spaced starting memories from 5.5 through 9.2.
3. For every start, record the full memory range over updates 6,145–8,192 and again over 24,577–32,768. Assign the low or high regime only when both ranges lie wholly on the same side of 7.5. Repeat with an independent complex-exponential implementation of the field map. Both windows gave the same classification for every run within each implementation. [Outcomes for every run](outcomes.json)

The binary64 runs agree on 3,019 of 3,088 dose cells. The two highlighted rows have eleven switches in the reference run. Keep just the cells where both runs agree and U04 shows seven, U05 one. The grid spacing is about 0.01927. The two late windows measure a finite-time regime. Finer bands and later transitions remain possible.

## Precision near the boundary

L01 and U04 had the most switches after excluding binary64 disagreements. For each, we repeated the 21 consecutive doses with zero-based indices 62 through 82, using the exact same binary64 starting coordinates, doses and coefficients. We increased the working precision five times. Every replay counts switches between neighboring samples on this interval.

| Field | 160 bits | 320 bits | 640 bits | 1,280 bits | 2,560 bits |
| --- | ---: | ---: | ---: | ---: | ---: |
| L01 | 9 | 7 | 5 | 5 | 5 |
| U04 | 7 | 5 | 7 | 7 | 7 |

All 42 labels matched at the three highest precisions. Eight of these outcomes differ from labels on which both binary64 runs agreed. One example is L01 at index 70, where both gave the upper regime and all three highest precisions gave the lower one. A separate 160-bit replay of **all 256 donor swaps** matched the donor's regime in every case. The replay keeps the midpoint after each update, so the agreement tests numerical repeatability. It supplies no interval enclosure of the exact orbit. [Labels, dose values and donor outcomes](precision.json)

## Reproduce

From the repository root, run the following with Node.js 18 or later, ImageMagick 7 and Python with NumPy on the path. Set the environment variable `PYTHON` to use a specific interpreter.

```sh
node research/art_20261122/memory-transplant.mjs
```

Both binary64 calculations start from the published donor coordinates. Node.js 18.16.0 supplied the reference labels and NumPy 2.5.5 supplied the replay. To reproduce the precision table and donor check, run `python research/art_20261122/memory_transplant_precision.py` with python-flint 0.9.0. [Reference calculation and figure](../../../research/art_20261122/memory-transplant.mjs) · [Independent binary64 implementation](../../../research/art_20261122/memory_transplant_replay.py) · [Higher-precision replay](../../../research/art_20261122/memory_transplant_precision.py)

The donors themselves depend on about 33,000 preceding chaotic updates. Repeating their sampling loop in Edge 158.0.4370.33 reproduced all 48 archived coordinates bit for bit. Node.js 18.16.0 gave different snapshots with the same arithmetic expressions. The optional `--check-donors` flag on the Node command compares a runtime's regenerated coordinates with the archive. [Map update and donor sampler](../../../research/art_20261122/memory_transplant_map.mjs)

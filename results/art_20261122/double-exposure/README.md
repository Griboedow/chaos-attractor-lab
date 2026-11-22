# One photograph, two histories

The two middle panels are identical, pixel for pixel. Run the top photograph’s particles fifteen steps backward and they draw a seven-lobed rosette. Run them fifteen steps forward and they draw an eye. The bottom photograph contains the same number of particles in every pixel. Its particles disperse in both directions.

![Two rows of three exposures. A rosette, a grainy square and an eye appear above three grainy squares. The two center squares are identical.](two-histories.png)

[Erase one direction of time](#erase-one-direction-of-time)

## How to hide both pictures

The clock is the integer-grid version of Arnold’s cat map with $`p=2^{31}-1`$ sites per axis. One tick multiplies each position by the matrix below, reducing both coordinates modulo $`p`$. This permutes a finite grid, so every individual particle is periodic. The fifteen steps in either direction use exact integer arithmetic.

```math
A=\begin{pmatrix}2&1\\1&1\end{pmatrix},
\qquad A^{-1}=\begin{pmatrix}1&-1\\-1&2\end{pmatrix}.
```

We selected the two endpoint images in advance. The past mask follows the two fold curves of the [archived rosette](../rosette/README.md#where-the-fourteen-cusps-come-from). The eye is a drawing made from two lids, an iris and a pupil. Both masks use 1,280 × 1,280 cells.

Sample independent, uniform present positions. Keep a position only if it lands inside the rosette at step −15 **and** inside the eye at +15. Repeat until 262,144 particles pass. For the second row, independently redraw each particle’s two integer coordinates within its original pixel on the 320 × 320 grid. Use the same dynamics, camera and brightness scale for both rows.

| Fraction inside the mask | Past rosette | Future eye | Both |
|---|---:|---:|---:|
| Selected particles | 100% | 100% | 100% |
| Subpixel positions redrawn | 11.23% | 9.12% | 1.04% |
| Mask area, or product of areas in the last column | 11.12% | 9.16% | 1.02% |

All 102,400 present pixel counts match exactly. At that resolution the selected cloud’s Pearson chi-square per degree of freedom against uniform density was 0.9993. An independently sampled uniform cloud gave 1.0073. Redrawing the subpixel positions creates a different physical state with the same photograph. At the two endpoints its mask-hit rates are close to the masks' area fractions.

Schulman studied two endpoint conditions and perturbations in the cat map, measuring how the coarse-grained entropy changes over time. Here the endpoint conditions specify two recognizable drawings and the control preserves every visible count at the middle time. [Schulman, *Causality Is an Effect, II* (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8229640/)

## Erase one direction of time

Start from the selected particles again. In the top row, give each particle a random displacement along the forward-stretching direction **inside its present pixel**. Fifteen steps backward still reveal the rosette. Fifteen steps forward give grain at the displayed resolution. Turn the displacement direction by 90 degrees and the eye survives instead. The two middle panels still match the original exposure in all 102,400 pixels.

![Two identical middle exposures. Moving along one direction preserves the past rosette alone. Turning that direction by ninety degrees preserves the future eye alone.](directional-eraser.png)

| Subpixel displacement | Inside the past rosette | Inside the future eye |
|---|---:|---:|
| Along the forward-stretching direction | 99.95% | 9.34% |
| Along the backward-stretching direction | 11.32% | 99.97% |
| Area of each mask | 11.12% | 9.16% |

Zero is an allowed displacement. In the top row 723 particles stayed exactly where they were, and in the bottom row 691 did. Those particles still hit both original masks. Among particles that actually moved, 9.09% hit the faded eye in the top row and 11.08% hit the faded rosette in the bottom row.

The two lattice directions are $`A^{10}(1,0)=(10946,6765)`$ and $`A^{-10}(0,1)=(-6765,10946)`$. They approximate the expanding directions of the map and its inverse. In each row, choose uniformly among all integer multiples of that row's vector that keep the particle in its pixel. Forward evolution magnifies a displacement along the first vector. Backward evolution magnifies a displacement along the second. This accounts for which drawing disappears. [Exact counts and vectors](directional-eraser.json)

## Reproduce

Run `python research/art_20261122/double_exposure.py` with NumPy, Pillow and Matplotlib. The script regenerates the six-panel plate, checks both forward/backward round trips for every point, and compares the two central histograms. [Source and mask definitions](../../../research/art_20261122/double_exposure.py) · [Counts, seeds and four grid resolutions](measurements.json)

Run `python research/art_20261122/directional_eraser.py` to regenerate the second plate from the same initial sample. The script checks each particle's present pixel before counting hits at the two endpoints. [Directional displacement code](../../../research/art_20261122/directional_eraser.py)

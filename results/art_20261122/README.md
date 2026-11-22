# Mathematical art from four dynamical systems

| Study | How to read the picture |
| --- | --- |
| [One rosette, two motions](rosette/README.md) | A rigid rotation and a chaotic map expose the same seven-lobed density. In the animation, moving each sample halfway towards its successor reveals which clock drew it. |
| [A nine-color Newton weave](newton/README.md) | One circuit of the torus wraps nine successive magnifications of a basin boundary. Five different poles send a common small displacement to five different roots. |
| [A mosaic of straight geodesics](weave/README.md) | The eight straight rays on a surface of genus fifteen divide it into 140 tiles. At equal wave weights sixteen more rays appear. One new line in the cosine plane accounts for the jump to 280 tiles. |
| [The seventeen-point crown](rings/README.md) | A single point paints seventeen lobes. The places where the map folds predict the sharp tips. Inside each tip three nearby preimages merge into one. |

Each study gives its equations, settings, evidence, source code and references. [Main atlas](../../README.md)

## Where to look next

**Five transparent maps of one Newton cycle.** Give all five poles the same displacement. Paint one transparent layer per pole wherever that start reaches the nine-cycle. Use angle horizontally and logarithmic distance vertically. Two layers already overlap. At radius 10⁻¹⁴ and angle 27.5°, direct interval iteration puts poles 2 and 5 in the cycle and the other three at distinct roots. Can three layers overlap anywhere within radius 1.79568 × 10⁻¹⁴? One certified patch would settle it. [Pole numbering and verifier](newton/README.md#five-poles-five-different-roots) · [Hubbard, Schleicher and Sutherland on Newton basins](https://www.math.stonybrook.edu/~scott/Papers/Newton-HSS.pdf)

**One orbit or seven?** In Case 180844, rotate the two archived periodic curves by one seventh of a turn and overlay them. With the best allowed time shifts, the longer curve matches itself to within 2 × 10⁻¹¹ coordinate units. For the shorter curve the maximum discrepancy is 0.130. Draw all seven rotations in distinct colors. Do they trace seven separate periodic orbits? Recompute the shorter branch near the jump in its recorded period to check whether the continuation switched branches. [Orbit data](../shape_search_20261121/flows/nested-skeletons.json) · [Wulff and Schebesch on following symmetric orbits](https://doi.org/10.1137/050637170)

**A flash on a surface with fifteen handles.** Send a wave in all tangent directions from one point. Color each patch by the route that reaches it first. A seam marks two equally short routes around the handles. When the straight-line mosaic jumps from 140 tiles to 280, does the seam network reconnect or slide smoothly through the change? Compare the same source at nearby weights, refine the surface mesh and trace the two competing geodesics back to the source. [Surface and mosaic](weave/README.md) · [Mitchell, Mount and Papadimitriou on shortest paths](https://www.cis.upenn.edu/~cis6100/discrete-geodesic-87-Mitchell-Mount.pdf)

**Two knots, one coloring rule.** Unroll the four-strand cores of Cases 56014 and 114425. Both knots have determinant one. Number the colors from 0 to 4. At a positive crossing the upper strand keeps color x and the lower strand changes from y to 2y − x modulo five. Reverse this update at a negative crossing. Of the 625 starting palettes, only the five solid colors close on the first knot. On the second, 20 more close and each uses all five colors. Follow one of those palettes to the seam, where the returning colors must match the starting ones. [Exact braid words](../shape_search_20261121/braids/certificates.json) · [Kauffman and Lopes on linear knot colorings](https://arxiv.org/abs/1708.01932)

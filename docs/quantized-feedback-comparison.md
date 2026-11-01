# Projection dynamics and quantized feedback

The finite-site RAAR examples must also be compared with digital filters,
sigma-delta modulators, overflow maps and quantized feedback systems. Switching
among discrete branches, irregular bounded motion and symbolic itineraries are
established features of this literature. They are not novelty claims for our
[circle examples](raar-constellations.md) or
[sphere survey](spherical-constellation-search.md).

## Sources inspected

| Source | Material inspected | Relevance and limit |
|---|---|---|
| [Wu and Chua, UCB/ERL M92/18 (1992)](https://digicoll.lib.berkeley.edu/record/137497/files/ERL-92-18.pdf) | Abstract, introduction and Section 2; PDF page 6 also checked visually | A second-order filter with overflow has three symbolic regions in the stated parameter range. This is a piecewise-affine map, not an identification of our projection map. |
| [Kocarev and Chua, UCB/ERL M92/74 (1992)](https://digicoll.lib.berkeley.edu/record/134463/files/ERL-92-74.pdf) | Abstract and Section 2 through the statement of Theorem 1; PDF pages 6–7 checked visually | Extends the symbolic description to other parameters. Equations 1–2 give the map and its linear part explicitly. We have not audited the complete proof of its chaos classification. |
| [Farrell and Feely, *Application of Non-Linear Analysis Techniques to Sigma-Delta Modulators* (1997)](https://mural.maynoothuniversity.ie/id/eprint/9754/1/RF-Application-1997.pdf) | Extracted text of PDF pages 1–4, including the geometric-analysis section | Discusses affine regions, trajectories on parabolic curves, switching and chaotic modulators. Figures and the poorly extracted band-pass equation were not reviewed visually; no exact comparison with that equation is claimed. |
| [Feely, *Nonlinear dynamics and chaos in Sigma-Delta modulation* (1994)](https://doi.org/10.1016/0016-0032(94)90092-2) | Publisher abstract and bibliographic information | Establishes an earlier review of nonlinear and chaotic converter dynamics; the full article was not inspected. |

For the overflow map displayed by Kocarev and Chua,

```math
G(x,y)=(y,-x+ay+2s(x,y)),
```

the integer-valued selector keeps the second coordinate in its prescribed
interval. On each open selector region the Jacobian is the same matrix
`[[0,1],[-1,a]]`. Wu and Chua use the corresponding three-region description
when `b=-1` and `|a|<2`. These observations concern the actual equations in the
inspected reports, rather than resemblance between pictures.

## A limited exact exclusion

There is a simple obstruction to an **ambient affine conjugacy on an open
set** between our RAAR maps and a piecewise-affine map. This excludes a
particular kind of coordinate disguise; it does not establish historical
novelty or topological inequivalence of their invariant sets.

Let `H(x)=Mx+c`, with nonsingular constant `M`. If
`H F H^{-1}=G`, then

```math
DG(Hx)=M\,DF(x)\,M^{-1},\qquad
\operatorname{tr}DG(Hx)=\operatorname{tr}DF(x).
```

On a fixed nearest-site branch of the circle map, let
`r=|(2b-x)_xy|`. Its trace is

```math
\operatorname{tr}DF(x)=3\beta-\frac{\beta}{r}.
```

For the unit-sphere map, with `r=|2b-x|`, it is

```math
\operatorname{tr}DF(x)=3\beta-\frac{2\beta}{r}.
```

For `0<beta<1`, either expression is nonconstant on every nonempty open
subset of a regular branch. A piecewise-affine map has constant Jacobian on
each open branch. Suppose its switching boundaries have empty interior and
form a locally finite closed union, as in ordinary finite quantizers and the
overflow maps above. Any proposed open-domain conjugacy contains a smaller
open set lying inside regular branches of both maps. Trace invariance then
gives a contradiction.

The statement applies to equal-dimensional affine models, including affine
state extensions of an overflow filter. A two-state filter and our three-state
map do not admit a nonsingular ambient affine change of coordinates in the
first place. The argument does **not** exclude a reduction restricted to a
lower-dimensional invariant set, a nonlinear conjugacy, a return-map relation,
a symbolic factor, or a nonlinear quantized-feedback model.

## What remains open

The closer comparison remains finite-alphabet phase and bit retrieval, whose
normalization branches are themselves nonlinear. The trace argument does not
exclude those algorithms. Their connection to our circle geometry is detailed
in the [existing comparison](raar-constellations.md#prior-art-and-equivalent-terms).

The present extension therefore supports two narrow conclusions: quantized
algorithmic chaos is established prior art, and ordinary piecewise-affine
filters are not affine copies of these RAAR maps on an open state domain.
No exhaustive classification of quantized-feedback models or nonlinear
equivalences has been completed. The [review record](../results/spherical-constellations/priority-review.json)
preserves source hashes, inspected pages and retrieval failures.

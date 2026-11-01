# Spatial RAAR: tangent-seed evolution and complete preimages

All six declared seed experiments, both numerical inverse searches, both interval-Sturm certificate revisions, independent nested-reflector certificates and all visual revisions are retained. Earlier RAAR packets are unchanged. Source paths in reports are provenance; use the local folders mapped in provenance.json.

The final figures are in `render` and `inverse-render`. Earlier versions retain their original annotations; the corrected final rendering scope in `render/analysis.json` supersedes earlier image-resampling wording. Every point contributes to the histograms. The final GIF uses native histogram pixels after a common empty-margin crop.

The near-cloud search found one numerical preimage per target. The subsequent global search and interval certificates found exactly two or four regular-domain preimages. Additional inputs are off the finite reference cloud: noninjectivity of the full map does not establish noninjectivity on an attracting set.

Run `python research/verify_raar_folding.py` from the repository root. See [the report](../../docs/raar-folding.md) and [the manuscript](../../papers/raar-folding/manuscript.pdf).

# Strict trapping for the spatial-circle RAAR instance R1051

The frozen full run processes **21,405,802** source cubes, accepting **23,055,541** dyadic leaves after **235,677** eight-way subdivisions. Every image, enlarged by the exact cube margin `2^-35`, is contained in the same occupied region. No new target cell is required. Both projection-axis arguments are uniformly bounded away from zero by more than `0.0118`.

The seven previously certified symbolic windows and the located period-7/49 root boxes are inside the region. The compactness/attraction and entropy arguments are in the [manuscript](../../papers/raar-trapping/manuscript.pdf) and [current report](../../docs/raar-trapping.md). Historical uniqueness, transitivity and equality with the finite visual cloud remain unproved.

## Contents

- `candidate/`: exact protocol and lossless gzip of the original packed mask. The sorted key file is reconstructed from the mask; both restored NPY files must match their original SHA-256 hashes. The uncompressed mask and keys occupy about 274 MB together.
- `runs/raar-trapping-audit-round2-part0..3`: the complete production calculation, including source snapshots, runtime declarations, all chunk records and all returned new-cell arrays. The four ranges partition the complete sorted key list.
- `runs/raar-trapping-closure-audit01`: separate mask/key, exact-window/cycle, trigonometric-bound and complete range accounting. It also recomputes the seven-cycle's negative determinant from its earlier interval monodromy.
- `runs/raar-trapping-audit-round1-part0..3`: original raster control flow. Three parts finish; the fourth is deliberately stopped after the complete optimized calculation passes. Its five-million-source checkpoint is retained and the unfinished in-flight chunk is excluded. This incomplete round is not counted as a second complete proof.
- `runs/raar-trapping-round-comparison01`: all **21,054,351** completed common sources agree exactly in accepted-leaf accounting, axis bounds and returned keys. Both control flows share their interval arithmetic and geometric enclosure.
- `runs/*validation*`: elementary arithmetic, high-precision geometric and exact-rational raster controls, including failed revisions. These selected checks are not substitutes for exhaustive interval execution.
- `runs/raar-trapping-render01`: a failed camera-dimension assertion before figure output, explicitly labeled incomplete. `render02` contains the reviewed final figures and their actual data. Both new figures were drawn on **2026-11-01**.
- `proposals/`: available source and protocol files from the numerical search, including failures. Large proposal-only arrays remain local and are explicitly listed, sized and hashed in `proposal-omissions.json`. All proof inputs are reconstructible without them.
- `priority.json`: bounded additional search and primary-source inspection, including Fannjiang and Zhang (2020). Third-party publication PDFs are not included.

The proposal uses sampled curvature to locate a candidate. Its numerical forward closure is not the proof. The strict run uses analytic geometric Taylor bounds and outward arithmetic on every source cell. The exact interval construction is described in the manuscript.

The outward kernels compile with cache disabled. Earlier development outputs retain their original sources and limitations. The final optimized kernel differs from the original only by stopping a target scan once subdivision is already necessary; it still scans all unresolved targets at maximum depth.

## Verification

Use the complete [reproduction commands](../../papers/raar-trapping/README.md). A fast verification restores inputs, checks every saved file and source range, audits the existing symbolic dependency and replays the finite figure data:

```powershell
python research/verify_raar_trapping.py --out outputs/trapping-verification
```

This does not independently regenerate all interval image enclosures. The four full batch commands in the reproduction guide do. The proof assumes the stated IEEE-754 basic-operation/square-root behavior and runtime preservation; it is not a formal proof-assistant certificate.

The historical symbolic packet is an explicit dependency identified by its manifest hash. Its frozen contents and earlier manuscript remain unchanged. `manifest.json` identifies every file in this new packet; `provenance.json` records copied source bytes. Current documentation and drawing captions remain outside the frozen numerical manifest.

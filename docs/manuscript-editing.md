# Manuscripts and historical publication records

The complete articles are edited directly in `papers/<study>/manuscript.md`.
They contain the proofs, equations, tables, references and links to the original
PNG figures. No document compilation is required. Display equations use GitHub's
`$$` math blocks; numerical results remain in the linked evidence packets.

Run from the repository root after editing:

```bash
python tools/check_manuscripts.py
python tools/render_dates.py sync
python tools/render_dates.py check
```

Use Markdown for new manuscripts and reports. Keep figures in the established
image formats and style. Do not commit PDF exports or maintain a second prose
source in Python or LaTeX. Numerical code and machine-readable proof records
retain their existing formats.

## Conversion and provenance

On 2026-11-01, all ten articles were transferred from their original authoring
sources to Markdown. The nine structured manuscripts retain the text and tables
actually consumed by their builders, including computed values. The motor article
retains its abstract, proofs, appendices, numbered equations and bibliography.
Only presentation and repository navigation changed.

The [conversion receipt](../provenance/markdown-migration.json) records the source
hashes, content inventories and hashes of the initial Markdown editions.
Subsequent edits need not match those initial hashes.

The existing `build.json`, `review.json`, `revision.json` and `manifest.json`
records inside paper directories describe the **historical publication editions**,
not the current Markdown files. They are retained unchanged. Frozen evidence
packet READMEs also retain their original bytes and may refer to archived
publication paths; use the current study reports or [manuscript index](../papers/README.md)
for reading links. Their page counts,
rendered-page hashes and PDF reviews must not be presented as reviews of a later
Markdown edit. The removed artifacts and their replacements are listed in the
[publication archive](../provenance/publication-archive.json), pinned to commit
`856ac3709801b0e42bb796e55c0e4fc2d8a11e04`. The concurrently published fold-membership
article uses the per-file source commit `2ff2973cfdfaa9bac23a7e7b4cf1ff22439de228`.

## Replaying a frozen numerical packet

Some original verifiers require the exact file inventory of a frozen packet,
including its former figure exports. Their checks and manifests remain unchanged.
Restore only the required historical presentation files before using such a
verifier, for example:

```bash
python tools/publication_archive.py --prefix results/raar-trapping --materialize
python research/verify_raar_trapping.py --out outputs/trapping-verification
```

Restoration is explicit. Files are checked against the original SHA-256 and read
from local Git objects, with a pinned GitHub download as a fallback. The default
archive cache is `D:/Codex/chaos-publication-archive`; `--output` selects another
external directory. Without `--materialize`, nothing is written into the working
tree. Restored publication files are ignored by Git. For omitted numerical arrays,
use the separate procedure in [data storage](data-storage.md).

To reproduce the original PDF build or publication audit, use a separate checkout
of the relevant pinned commit on D:. That revision contains the original authoring sources,
builders, auditors and publication file inventory. These historical audits are
not checks for the current Markdown workflow.

# Public repository rendering check

The separate [Chaos Attractor Lab repository](https://github.com/Griboedow/chaos-attractor-lab) was checked on 30 October 2026 at immutable initial commit [`a514509`](https://github.com/Griboedow/chaos-attractor-lab/commit/a5145092d86ff5ad51d42bc5e1a7a4eb64d2d802). GitHub reports **public** visibility and **MIT** licensing. The published license text matches that commit's file and retains Copyright (c) 2026 Griboedow.

This is a new-repository check. The older [formatting history](report-format-audit.md) records checks of the Chua repository at its stated revisions; those successes are not counted again here.

## Actual public-page checks

The pages were loaded anonymously in installed Microsoft Edge through headless Playwright, with JavaScript and MathJax enabled, at a 1320 × 1000 desktop viewport.

| Page | Rendered math elements | Result |
| --- | ---: | --- |
| [README](../README.md) | 0 | Pass |
| [Research findings](research-findings.md) | 0 | Pass |
| [Synchronematic research](synchronematic-research.md) | 8 | Pass |
| [Coupled optical operators](coupled-operators-research.md) | 0 | Pass |
| [Optical lifetimes to 100,000](optical-lifetime-extension.md) | 0 | Pass |
| [Optical periodic-cycle continuation](optical-cycle-continuation.md) | 15 | Pass |
| [Thermal-flutter source audit](thermal-flutter-source-audit.md) | 16 | Pass |
| [Optical-order constraint](coupled-operators-order-audit.md) | 36 | Pass |

All **75 math elements across eight pages** rendered without unresolved formulas, error/red command tokens or math-container overflow. A separate scan outside math and code found no raw dollar delimiters or TeX remnants. Zero math elements means the page uses prose, tables and code-formatted labels; it does not mean its content was skipped. No formatting repair was needed.

Every visible article image loaded. The complete optical-lifetime figure, findings overview, logistic feedback law and six aligned operator equations were visually inspected. Two browser screenshots of the synchronematic GIF, separated by 1.1 seconds, show different camera orientations: animation was checked, not just the existence of an image URL. Hidden GIF fallback images were recorded separately. The README's **Read the findings** link was clicked and opened the matching report at the same immutable revision. These supporting image and navigation visits are not extra full-page checks.

## Relative paths and limits

All **372 relative Markdown path occurrences across the commit's 63 Markdown files** resolve in its Git tree. The matcher excludes fenced code/math and inline code before interpreting links. This broader path check does not imply browser inspection of all 63 pages. External references, general fragment anchors, other browsers, mobile layouts and later edits remain outside scope.

These checks concern presentation and public availability. They do not establish scientific novelty, validate every numerical result or replace a separate fresh-clone runtime check. This note was prepared after the inspected initial commit and is not one of its eight checked pages.

The separate [fresh-clone runtime check](../provenance/public-clone-validation.json) subsequently passed at the same initial commit. It remains a distinct numerical packaging record and is not included in these browser-check counts.

The [machine-readable record](../provenance/github-rendering.json) retains exact page URLs, image diagnostics, repository/license evidence, GIF frame hashes and the clicked-link result. Browser UI screenshots remain in ignored local outputs and are not redistributed under the repository's MIT license.

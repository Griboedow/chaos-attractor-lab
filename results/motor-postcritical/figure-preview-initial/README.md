# Superseded visualization preview

These PNG files reproduce the first visual-QA rendering exactly. Their hashes match the historical `../figure-manifest.json`. The final public figures are identified by `../figure-manifest-v2.json` and live under `docs/assets/`.

The first heatmap used interpolation across discrete pair rows, and the final-window trace joined a sparse trajectory sample. Visual QA replaced these with a noninterpolated pair heatmap and smooth phase-aligned periodic-orbit traces. The numerical data and conclusions did not change. The initial render also encountered a Windows console encoding error only when printing an absolute path after all images and its manifest had been saved; the current renderer prints repository-relative paths.

# Unexpected-search checkpoint, 2026-11-01

This packet supports [the tetrahedral flow](../../docs/tetrahedral-flow.md) and [the scalar acceleration-memory investigation](../../docs/anderson-memory.md). Both are numerical chaotic-attractor candidates; neither is assigned a historical-novelty probability.

The packet preserves two complementary approaches, including their failed candidates:

- `tetra-v1`: all 1,024 parameters, starts, growth windows and scout clouds; 127 initial positives.
- `tetra-long`: the declared 36-row selection and longer smaller-step follow-up; 26 survivors.
- `tetra-spectra`: four candidates at three step sizes, full spectra and disjoint windows.
- `tetra-dop853`: independently implemented adaptive solver and Jacobian for two candidates.
- `tetra-basins`: fresh volume, local and stable-origin starts.
- `tetra-render`: 40,000-unit dynamics, 100,000-point samples, camera previews and twelve-rotation support diagnostics. The two dense `.npy` arrays remain scratch artifacts and can be regenerated.
- `anderson-v1`: the negative 4,096-row contraction screen.
- `anderson-v2`: a separate 4,096-row extended screen and all 132 initial positives.
- `anderson-long`: all 132 candidates observed for 500,000 steps; failures retained.
- `anderson-spectra`: full spectra and independent exact linear-reduction controls.
- `anderson-precision`: 80/110-digit calculations.
- `anderson-basins`: interval verification of unique scalar roots and 626 fresh starts.
- `comparison`: 108 known regimes plus six previous repository candidates; raw and covariance-whitened descriptors and local alignments.
- `visuals`: images, rigid rotation, histograms, sample clouds and exact rendering provenance.
- `source-original`: the exact calculation sources and their reused repository dependencies.

`facts.json` supplies machine-readable accounting and parameters. `environment.json` records the environment and corrected implementation failures. `reference-inputs.json` pins external inputs already stored elsewhere in this repository. `manifest.json` covers this packet, the portable scripts, the two reports and their published assets, excluding the manifest itself. Shared landing-page text is intentionally outside the manifest.

All new scratch and repository work was on drive D. See [portable reproduction instructions](../../research/unexpected_search/README.md). The archived finite calculations do not prove infinite-time chaos, topological novelty or that nobody has previously seen these geometries.

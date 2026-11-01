# Reference clouds and bounded novelty audit

[Report and scope](../../docs/attractor-novelty-audit.md) · [New vortex image](../../docs/gaussian-vortex.md)

`reference/trajectories/` contains 108 complete finite point clouds. Counts are 89 eligible dysts defaults, four independently transcribed public-atlas equations, four published 3D Hénon regimes, nine successful additional Halvorsen parameters, and two kicked-top tails imported from an earlier frozen repository study. Failed Halvorsen attempts and the excluded fixed-point kicked-top outcome are explicit in their records.

`reference/protocol.json` identifies every included/excluded catalog model and its defaults. The `*-runs.json` files give the remaining protocols. Upstream source and licenses are retained for attribution; public README novelty claims are not this study's conclusions.

`comparison/summary.json` contains the rotation/translation/reflection/uniform-scale screen, five nearest fits per candidate, within-trajectory controls and subsample rankings. `affine-summary.json` adds covariance whitening to remove nonsingular affine coordinate changes. `aligned-clouds.npz` preserves the fitted comparison samples; `plotted-clouds.npz` and `render.json` reproduce the separate figure.

`literature-review.json` records searched terms, inspected primary sources and outstanding limitations. Existing algorithmic and motor studies have separate priority audits. These geometric scores do not prove topological or historical novelty, and local ICP fits are not certified global optima.

Both packets share [the Gaussian-vortex manifest](../gaussian-vortex/manifest.json) and `research/verify_gaussian_packet.py`. Original executed reference-build scripts are in `executed-source/`; portable scripts in `research/` accept an output directory. All numerical comparisons can be repeated from the committed clouds without accessing the internet.

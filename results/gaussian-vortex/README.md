# Gaussian-vortex evidence packet

[Illustrated report](../../docs/gaussian-vortex.md) · [Novelty audit](../../docs/attractor-novelty-audit.md)

`candidates.json` defines the equation parameters and their order. `provenance.json` records the environment and selection protocol. `search-summary.npz` retains every one of 1,024 scout rows; `followup-summary.npz` retains all 45 follow-ups, including four failures of the stated growth gate. Their galleries remain available. Full exploratory tails are indexed in `local-large-files.json` rather than committed.

`long-runs.npz` contains six full-spectrum histories, ordered by candidate 77/233/738 and step 0.0125/0.00625. Every history discards 2,000 time units and observes 12,000; windows are 500 units. Each retained tail has 50,000 states, sampled every four steps at the end of the run.

`render-run.npz` records the separate 30,000-unit primary trajectory, its thirty full-spectrum windows and final 100,000 states. `render-histogram.npz` counts all 4.8 million displayed states. `render.json` records the raw trajectory hash and exact camera/color choices. The large full cloud remains locally indexed and is reproducible by the portable renderer.

`fresh-starts.npz` and `neighbors.npz` retain every declared initial state and parameter row, full spectra and tails. `controls.json` records independent algebra, exact linear solutions, local integration errors and adaptive full spectra. `regular-controls.json` continues all four regular fresh-start outcomes; they show numerical period 63.83832949, not a validated periodic-orbit theorem. `summary.json` is derived from these arrays.

`executed-source/` preserves the original search and rendering implementations, including their local path bindings. Portable entry points live under `research/`. `replay-checks.json` records their separate replay. `manifest.json` freezes this packet, the related `attractor-novelty` packet, dedicated reports, code, images, and the earlier input data used for comparison. Shared repository indexes are intentionally outside the manifest.

Run `python research/verify_gaussian_packet.py` from the repository root. This checks integrity and numerical consistency; it does not certify chaos or novelty.

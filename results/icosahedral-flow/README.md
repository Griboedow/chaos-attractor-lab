# Icosahedral-flow search checkpoint

All 256 scout rows, the five-candidate refinement (including failed case 85),
three long trajectories, two independent adaptive integrations, fresh basin
starts, actual-state rendering data, finite chirality diagnostics and geometric
comparisons are retained. Reiter's 1999 Figure 1 is independently replayed from
its published formula and passes all four rounded iteration checkpoints.

The 70-decimal outward interval bounds certify absorbing balls for the stored
coefficients. Other chaos, attraction and symmetry claims are numerical.
Historical uniqueness remains unresolved. See ../../docs/icosahedral-flow.md.

Dense 1,280,000-state arrays are regenerated using icosa_long.py. Their strides,
byte hashes, 80,000-state samples and unmodified density histograms are archived.
The GIFs rotate a fixed subset of one trajectory; no group copies are added.

Portable scripts: ../../research/icosahedral_search/README.md.
Manifest paths are relative to the repository root; the manifest does not hash
itself. Existing earlier packets are treated as immutable inputs.

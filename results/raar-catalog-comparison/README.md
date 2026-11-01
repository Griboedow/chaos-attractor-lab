# Bounded RAAR geometry comparison

The existing million-state R1051 sample is compared with the 108 frozen regimes in results/attractor-novelty/reference. No reference system is reintegrated. Both shape-preserving and covariance-whitened representations retain full rankings, four resampling controls, temporal-half comparisons and three local alignments. The plotted points are actual stored states. The initial rendering and corrected caption are both retained.

This is a finite empirical comparison, not a proof of dynamical inequivalence, historical novelty or chaotic attraction. Local alignment does not certify an optimal residual. Run research/verify_raar_catalog.py. Reproduce with research/compare_raar_catalog.py --out NEW_DIRECTORY, then research/render_raar_catalog.py --source NEW_DIRECTORY --out NEW_RENDER_DIRECTORY.

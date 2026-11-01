# Paired sphere / finite-site RAAR survey

See the [scientific report](../../docs/spherical-constellation-search.md) and
[quantized-feedback comparison](../../docs/quantized-feedback-comparison.md).
The archive retains six completed run directories: the paired screen, all-seven
follow-up, two long-run batches covering every flag, dense rendering of two
cases, and the exact rational period-two checks.

`inventory.json` lists all 59 original run files. Thirteen arrays totaling
93,968,695 bytes are optional cache payloads. Their exact inputs and compact
metrics are retained under distinct names in `results/compact-data`; see
`omissions.json` and the [storage instructions](../../docs/data-storage.md).
No large payload was committed and later silently replaced by a thin trajectory.

## Reproduction

Use a new external output directory for every calculation. The original circle
survey supplies the exact same parameters, finite sets and identifiers:

```text
python research/search_raar_constellations.py --out D:/Codex/new-circle-scout --power 9
python research/search_spherical_constellations.py --source D:/Codex/new-circle-scout/screen.npz --out D:/Codex/new-sphere-scout
python research/follow_spherical_constellations.py --source D:/Codex/new-sphere-scout/screen.npz --out D:/Codex/new-sphere-follow
python research/extend_spherical_constellations.py --source D:/Codex/new-sphere-scout/screen.npz --out D:/Codex/new-sphere-long01 --ids 1790 1119 --selection "Structured spatial candidate and volumetric comparison"
python research/extend_spherical_constellations.py --source D:/Codex/new-sphere-scout/screen.npz --out D:/Codex/new-sphere-long02 --ids 668 1086 1346 1694 1788 --selection "Remaining five scout flags"
python research/render_spherical_constellations.py --source D:/Codex/new-sphere-follow --out D:/Codex/new-sphere-render --images docs/assets/experimental/new-sphere-replay --ids 1790 1119
python research/check_spherical_cycle.py --source D:/Codex/new-sphere-render/case-1790-states.npz --site 6 --out D:/Codex/new-sphere-cycle
```

The first long batch predates the selection-description command-line argument;
its exact executed source is retained. Both batches call the same full-spectrum
implementation. A reproduced chaotic trajectory can differ between numerical
environments; historical byte hashes refer to the recorded environment and
original arrays, not a promise of platform-independent floating-point motion.

The paired screen input SHA-256 is
`da0a0c2a665a02c36bc408e931a5a771fa89447222bf9e3819f15b0d7c2eeb54`.
After restoring the sphere and circle optional inputs into the verified cache,
the replay/accounting check is:

```text
python research/verify_spherical_campaign.py --out D:/Codex/new-sphere-verification
```

The exact certificate uses rational arithmetic on the binary64 input fractions.
It proves the existence, uniform site admissibility and normal attraction of a
circle of two-cycles for Case 1790. It does not prove that the exact initial
condition enters its basin. The impressive finite portrait is retained as a
rejected transient. Historical novelty is unestablished for every case.

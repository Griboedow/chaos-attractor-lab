# Compact numerical records

These archives retain inputs and small numerical summaries from the optional
large payloads listed in [`provenance/generated-data.json`](../../provenance/generated-data.json).
A filename is the SHA-256 of its original full archive. Identical full archives
share one compact record.

Each retained array preserves its dtype, dimensions and bytes. The extraction
retains arrays up to 128 KiB and explicitly identified parameters, initial
states, seeds, field coefficients, cameras, step sizes, flags, statuses and
summary-rate windows. The registry enumerates every retained and omitted key;
the compact file has its own independent byte count and checksum.

These are partial archives, not replacements for the original files. A
calculation requiring trajectories or dense grids must recompute them or use
the [optional data cache](../../docs/data-storage.md). Frozen manifests continue
to refer to full historical payloads.

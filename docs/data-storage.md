# Numerical data storage

The working tree retains source code, equations, parameters, initial conditions,
compact numerical summaries, certificates and published figures. Large computed
trajectories, point clouds, histograms and intermediate arrays are optional cache
data. They are excluded from Git rather than committed after each calculation.

The [payload registry](../provenance/generated-data.json) records each omitted
file's original path, SHA-256, Git object, byte count, array shapes and retained
summary keys. The source snapshot is pinned to commit
`7468bb17bcc19db23103f3204079d0269cbd45e4`. [Compact arrays](../results/compact-data/README.md)
preserve small values and calculation inputs bit for bit; they are partial
archives and cannot replace a full trajectory in a numerical test.

Existing packet descriptions and frozen manifests describe the original full
calculation. Their recorded hashes remain unchanged. Reading the reports and
viewing the portraits requires no numerical download.

## Obtain a small checkout

```text
git clone --depth 1 https://github.com/Griboedow/chaos-attractor-lab.git
```

Removing files from the current tree does not remove them from historical Git
objects. A full-history clone still downloads the earlier large payloads. This
cleanup does not rewrite shared history.

## Recompute a study

Use the selected study's [reproduction instructions](reproduction.md) and its
recorded environment. Write new calculation outputs outside the checkout, on
D: where available. Parameters and initial states from an omitted archive are
listed in its registry entry and, when present, its `compact` file. Full source
snapshots and scripts remain in the original research packets.

A new numerical run need not reproduce the old binary bytes: chaotic
trajectories can separate under compiler or library changes. Use the original
recorded payload when a historical hash check is required.

## Recover a selected historical payload

The helper uses a local Git object when available, otherwise downloads only
the selected file from the pinned historical commit. It verifies byte count
and SHA-256 before accepting it. Nothing downloads implicitly during rendering.

```powershell
$env:CHAOS_DATA_CACHE='D:/Codex/chaos-attractor-cache'
python tools/data_cache.py list --group icosahedral-flow
python tools/data_cache.py fetch results/icosahedral-flow/icosa-long/sample.npz
```

On Windows the default cache is on D: when that drive exists. Elsewhere it is
under the user's cache directory. `--cache` selects a destination for a fetch;
set `CHAOS_DATA_CACHE` to the same directory when using the renderers.

To render the current presentation, fetch its required inputs only:

```text
python tools/data_cache.py fetch --portraits
python tools/render_portraits.py --output docs/assets/portraits --record-dates
python tools/animate_portraits.py --output docs/assets/portraits
python tools/check_portraits.py --with-data
```

Old calculation and manifest-verification scripts expect their original paths.
Use `--materialize` to place verified copies at those ignored paths:

```text
python tools/data_cache.py fetch --group icosahedral-flow --materialize
```

Some checks read other packets as well. Fetch the specifically requested input
or group. `fetch --all --materialize` restores the complete omitted numerical
snapshot and needs about four gigabytes; it is unnecessary for normal browsing.

## Keep subsequent commits small

Keep generated numerical payloads of 256 KiB or more outside Git. Exceptions
must identify an essential input or explicitly linked evidence. Preserve small
inputs, final metrics, provenance, code and the figures used in reports. Do not
replace a full array with a thinned array under the same filename or checksum.

Run `python tools/check_data_storage.py` before publishing. The normal
`python tools/check_portraits.py` checks published assets and recorded input
hashes without requiring cached arrays; `--with-data` also recomputes the
numerical rendering checks.

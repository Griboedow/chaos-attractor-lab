# Portrait conventions

The main README contains all current portraits from `docs/assets/portraits`.
Selected findings come first, followed by further portraits grouped by system
and the research index. There is no separate gallery page.
The [catalog](../records/portrait-catalog.json) identifies each numerical
source, coordinate selection and report. [Render records](../docs/assets/portraits/render.json)
give input hashes, sample counts, cameras and display operations. Archived
figures and numerical evidence retain their original bytes.

The [folder-change record](../records/layout-change.json) records the rename of
the record directory. Original render-source hashes remain unchanged in the
figure receipts. The checker accepts those sources only when reversing the
recorded path substitutions reproduces their exact original hashes; any other
source edit requires a new render receipt.

Static portraits are 1920 by 1080 pixels, with a white background, one blue density
palette, the same typography and fixed label positions. Panel labels identify
cases or coordinates; the scientific claim belongs in the surrounding text.
`Case 177`, for example, identifies a saved parameter set in the linked report;
it is not a measured quantity. Physical parameter labels include their meaning.

Existing density arrays are rebinned by summing square blocks. Only empty
outer margins are removed. State arrays use declared physical coordinates or
observables, a rigid orthographic camera and equal displayed units. Explicit
physical unit changes, such as `1000 I1`, are stated beside the figure. Paired
motor comparisons use a common spatial scale where declared.

Density uses `min(log1p(count)/q0.9995, 1)^0.72` with separate saturation per
panel. Brightness cannot be compared quantitatively between studies. Saved
adaptive knots are not uniform-time samples. Where continuous-flow points are
joined, straight segments are visual guides between consecutive saved states;
no additional trajectory states are inferred. Discrete map images are not
joined by curves. No symmetry replication, smoothing or synthetic image
generation is used.

Animated main-page portraits are 800 by 480 pixel lossless WebP files
in the same palette. A rigid camera rotates around the displayed vertical axis
through 600 views, each shown for 40 milliseconds: 25 frames per second and
one revolution in 24 seconds. Consecutive camera angles differ by 0.6 degrees. The
trajectory, labels, spatial bounds and intensity scale remain fixed for the
whole revolution. There is no simulated time evolution. The corresponding
full-resolution static portrait is embedded directly below each animation.
Every cataloged study has exactly one static portrait on the main page, with
its drawing date, coordinate description and evidence link.

[Animation records](../docs/assets/portraits/animation.json) give the source
hashes, retained states, orthogonal transforms, frame counts and durations.
The icosahedral inputs retain 320,000 states per case; Anderson uses its archived
100,000-state sample. RAAR, RK45 and smooth Yogi use the last 250,000 consecutive
states of their archived animation inputs. Finite-point RAAR uses the last 250,000
consecutive states of each two-million-state input, with the static PCA camera. Motor projections use the same
saved states as the static portraits. The [icosahedral input record](../results/portrait-presentation/README.md)
identifies the denser sample and its check against the frozen evidence.

```powershell
$env:MPLCONFIGDIR='D:/Codex/matplotlib-cache'
$env:PYTHONDONTWRITEBYTECODE='1'
python tools/render_portraits.py --output D:/Codex/portrait-review
```

For publication into the repository:

```text
python tools/data_cache.py fetch --portraits
python tools/render_portraits.py --output docs/assets/portraits --record-dates
python tools/animate_portraits.py --output docs/assets/portraits
python tools/render_dates.py sync
python tools/render_dates.py check
python tools/check_portraits.py --with-data
```

Use the NumPy and Matplotlib versions recorded in
[`results/icosahedral-returns/environment.json`](../results/icosahedral-returns/environment.json).
WebP encoding uses Pillow 13.0.0 with lossless compression.
The optional `--ids` argument selects named catalog entries for a preview.
Run the full catalog before publication so the final render record is complete.

Large input arrays are optional [cache data](../docs/data-storage.md). The
renderers read the external cache directly; they do not require copies in the
checkout. The normal checker, without `--with-data`, verifies published figures,
dates and input receipts without loading trajectories.

Dates accompany attractor portraits. Spectra, parameter scans, return maps and
other diagnostic plots do not receive drawing-date captions. A new rendering
links to the earlier registered portrait of the same study: the first recorded
date is preserved and a later rendering date is shown separately. These are
image records, not discovery dates. See [the date workflow](README-render-dates.md).

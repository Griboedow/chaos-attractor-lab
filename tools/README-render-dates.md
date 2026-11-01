# Drawing dates

The [dated gallery](../docs/attractor-gallery.md) and README display a date beside
every embedded local image. The machine-readable registry is
[`provenance/render-dates.json`](../provenance/render-dates.json).

For a new image, register the date immediately after the successful render:

```text
python tools/render_dates.py record docs/assets/new-candidate.png
```

This records today's date, the local timestamp including UTC offset, and the
image hash. If registering an earlier render, supply its actual date and source:

```text
python tools/render_dates.py record docs/assets/new-candidate.png --date 2026-11-01 --evidence results/new-candidate/render-log.json
```

Add the image and its evidence/status link to `docs/attractor-gallery.md`, and
optionally to the README. Then synchronize and check the captions:

```text
python tools/render_dates.py sync
python tools/render_dates.py check
```

The commands never modify image bytes. Re-registering an existing image preserves
its first date. Do not reset it for a new camera or higher resolution. The check
fails if a newly embedded image lacks a record or a caption is out of date.
The [repository instructions](../AGENTS.md) make these steps part of future work.

The initial backfill inspected complete Git history, including the original
`chua-gpu-atlas` repository for transferred figures. Each old entry includes its
earliest found image commit, author timestamp with timezone, path and blob ID.
These entries are explicitly labeled **First recorded drawing**: a commit proves
the picture existed then, but does not recover the exact time it was rendered.
The dates are not paper publication dates or historical discovery claims.

Frozen reports and evidence packets retain their original bytes and hashes.
Captions live in the current overview and dated gallery; their update does not
require redrawing a scientific figure or rewriting an archived manifest.

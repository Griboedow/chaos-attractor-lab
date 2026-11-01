# Portrait dates

Dates accompany the principal attractor portraits in the README, gallery and
current reports. Spectra, parameter scans, return diagrams and other diagnostic
figures do not receive drawing-date captions. All records remain in
[`provenance/render-dates.json`](../provenance/render-dates.json).

Register a new portrait immediately after its first successful render:

```text
python tools/render_dates.py record docs/assets/new-candidate.png
```

A new presentation of an existing study links to its earlier portrait:

```text
python tools/render_dates.py record docs/assets/portraits/new-view.png --source-portrait docs/assets/earlier-view.png
```

The record stores the actual rendering date, timestamp and image hash. The
caption preserves the earlier study portrait date and shows the current
rendering date separately when different. Re-registering the same asset does
not reset its first date. An earlier date requires `--date YYYY-MM-DD --evidence
<source>`; do not infer it from file-copy timestamps.

For an optional diagnostic record use `--role diagnostic`. Existing diagnostic
records are retained as provenance but excluded from captions.

```text
python tools/render_dates.py sync
python tools/render_dates.py check
```

These commands synchronize registered portrait captions and remove stale date
captions from diagnostic figures. New current `docs/*.md` pages are discovered
automatically. New main-page or gallery images must be registered. Archived
reports listed under `frozen_document_exclusions` keep their original bytes.

Historical records were reconstructed from complete Git history, including the
source `chua-gpu-atlas` repository. **First recorded portrait** is the earliest
found image-bearing commit date, with its original timezone; it is not an exact
recovered rendering time or a discovery date. The registry retains the commit,
path, timestamp and blob evidence.

See [portrait rendering](README-portraits.md) and [repository conventions](../AGENTS.md).

# Repository conventions

- Write repository content in English. Use ISO dates (`YYYY-MM-DD`).
- Keep generated trajectories, clouds, grids and intermediate arrays of 256 KiB
  or more outside Git, normally in `D:/Codex/chaos-attractor-cache` or a separate
  D: run directory. Retain code, parameters, initial conditions, compact numerical
  summaries, certificates and published figures. Do not replace full data with
  thinned data under an existing filename. Follow `docs/data-storage.md` and run
  `python tools/check_data_storage.py` before publishing. Large essential inputs
  require an explicit reason in the data registry. Do not re-add ignored caches.
- Every new attractor portrait or gallery entry must show its drawing date next
  to the image. Immediately after its first successful render, run
  `python tools/render_dates.py record <repository-relative-image-path>`.
  This records the actual local date and timestamp; use `--date` and `--evidence`
  when recording an earlier render with evidence. Never infer a drawing date
  from a clone, checkout, copy timestamp or scientific publication date.
- Dates belong to attractor portraits, not diagnostic plots. Register spectra,
  parameter scans and return diagrams with `--role diagnostic` when a record is
  useful; do not add date captions beside them.
- Keep the main README focused on key findings, their evidence level and links
  to child reports. Use concise academic prose; omit promotional language,
  conversational progress logs and repeated qualifications.
- Group main-page findings under descriptive section headings, with a
  subheading for each study. Use slow, smooth camera rotations (24 seconds per
  revolution at 25 frames per second) from `tools/animate_portraits.py` on the main page and link each
  to its static portrait. Keep trajectories and spatial scales fixed during
  rotation; retain static images in the gallery.
- Label stored parameter-set identifiers as `Case <number>` and explain their
  meaning beside the figures. Name physical parameters explicitly; do not put
  unexplained numbers above portraits.
- Use `provenance/portrait-catalog.json` and `tools/render_portraits.py` for the
  current README and gallery portraits. Keep the common layout and palette.
  Render from numerical data, identify projections and physical unit changes,
  and preserve archived originals. Use `--source-portrait` when a new
  presentation asset depicts an already registered study.
- Add new visual candidates to `docs/attractor-gallery.md`, including their
  evidence/status link. Keep paired cases identifiable. Preserve the first
  drawing date when changing the camera, resolution or rendering the same asset
  again. A genuinely new image asset gets its own record.
- When presenting a found attractor image in chat, include its drawing date
  from the same registry, in the user's conversation language.
- Before publishing gallery changes, run `python tools/render_dates.py sync`
  and `python tools/render_dates.py check`. Captions and the registry at
  `provenance/render-dates.json` must agree. Old dates reconstructed from Git
  remain labeled **First recorded drawing**, not exact rendering timestamps.
- Keep frozen numerical packets, images and their manifests unchanged when
  adding date captions. Use D: for temporary work and caches where possible.

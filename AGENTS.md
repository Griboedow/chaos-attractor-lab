# Repository conventions

- Write repository content in English. Use ISO dates (`YYYY-MM-DD`).
- Keep generated trajectories, clouds, grids and intermediate arrays of 256 KiB
  or more outside Git, normally in `D:/Codex/chaos-attractor-cache` or a separate
  D: run directory. Retain code, parameters, initial conditions, compact numerical
  summaries, certificates and published figures. Do not replace full data with
  thinned data under an existing filename. Follow `docs/data-storage.md` and run
  `python tools/check_data_storage.py` before publishing. Large essential inputs
  require an explicit reason in the data registry. Do not re-add ignored caches.
- Every new attractor portrait or main-page entry must show its drawing date next
  to the image. Immediately after its first successful render, run
  `python tools/render_dates.py record <repository-relative-image-path>`.
  This records the actual local date and timestamp; use `--date` and `--evidence`
  when recording an earlier render with evidence. Never infer a drawing date
  from a clone, checkout, copy timestamp or scientific publication date.
- Dates belong to attractor portraits, not diagnostic plots. Register spectra,
  parameter scans and return diagrams with `--role diagnostic` when a record is
  useful; do not add date captions beside them.
- Keep key findings first in the main README, followed by further portraits and
  the research index. State evidence levels and link to child reports.
  Use concise academic prose; omit promotional language,
  conversational progress logs and repeated qualifications.
- Keep the research index near the end of the main README, after the selected
  findings. `docs/research-index.md` is a compatibility link to that section,
  not a second copy of the index.
- Group main-page findings under descriptive section headings, with a
  subheading for each study. Use slow, smooth camera rotations (24 seconds per
  revolution at 25 frames per second) from `tools/animate_portraits.py` on the main
  page and embed the corresponding static portrait directly below each animation.
  Keep trajectories and spatial scales fixed during rotation.
- Label stored parameter-set identifiers as `Case <number>` and explain their
  meaning beside the figures. Name physical parameters explicitly; do not put
  unexplained numbers above portraits.
- Use `records/portrait-catalog.json` and `tools/render_portraits.py` for the
  current README portraits. Keep the common layout and palette.
  Render from numerical data, identify projections and physical unit changes,
  and preserve archived originals. Use `--source-portrait` when a new
  presentation asset depicts an already registered study.
- Add new visual candidates to the main README, including their evidence/status
  link. Keep every cataloged static portrait on that page exactly once; do not
  create a separate gallery page. Keep paired cases identifiable. Preserve the first
  drawing date when changing the camera, resolution or rendering the same asset
  again. A genuinely new image asset gets its own record.
- When presenting a found attractor image in chat, include its drawing date
  from the same registry, in the user's conversation language.
- Before publishing portrait changes, run `python tools/render_dates.py sync`
  and `python tools/render_dates.py check`. Captions and the registry at
  `records/render-dates.json` must agree. Old dates reconstructed from Git
  remain labeled **First recorded drawing**, not exact rendering timestamps.
- Keep frozen numerical packets, images and their manifests unchanged when
  adding date captions. Use D: for temporary work and caches where possible.

- Edit manuscripts directly in `docs/<study>/manuscript.md`; use Markdown for
  new manuscripts and reports. Do not commit PDF exports or duplicate prose
  sources in Python/LaTeX. Preserve historical review and manifest records with
  their original scope. Follow `docs/manuscript-editing.md` and run
  `python tools/check_manuscripts.py` before publishing manuscript changes.
- Use fenced `math` blocks for display equations so Markdown preserves TeX
  escapes, braces and subscripts. Protect inline TeX with GitHub's dollar/backtick
  syntax. Use upright math names for operators supported by the GitHub renderer.
  After publishing Markdown changes, open the
  changed pages on GitHub in a browser, wait for math rendering, and inspect
  formulas, tables, images and navigation. Check the published revision, not
  only a local preview. Fix rendering failures and repeat the browser check.

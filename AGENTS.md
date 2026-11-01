# Repository conventions

- Write repository content in English. Use ISO dates (`YYYY-MM-DD`).
- Every new attractor portrait or gallery entry must show its drawing date next
  to the image. Immediately after its first successful render, run
  `python tools/render_dates.py record <repository-relative-image-path>`.
  This records the actual local date and timestamp; use `--date` and `--evidence`
  when recording an earlier render with evidence. Never infer a drawing date
  from a clone, checkout, copy timestamp or scientific publication date.
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

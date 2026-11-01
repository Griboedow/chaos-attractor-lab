"""Check the editable manuscript format and local publication links."""
from __future__ import annotations

import json
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r"\]\(([^\s)]+)\)")


def main() -> None:
    failures = []
    manuscripts = sorted((ROOT / "papers").glob("*/manuscript.md"))
    expected = {v["paper"] for v in json.loads((ROOT / "provenance/markdown-migration.json").read_text(encoding="utf-8"))["manuscripts"]}
    actual = {p.parent.name for p in manuscripts}
    if expected - actual:
        failures.append(f"Missing migrated manuscripts: {sorted(expected - actual)}")
    # New manuscripts may be added without changing the historical conversion receipt.
    for path in manuscripts:
        text = path.read_text(encoding="utf-8")
        name = path.relative_to(ROOT)
        if not text.startswith("# ") or len(re.findall(r"^# ", text, re.M)) != 1:
            failures.append(f"{name}: expected one document title")
        if "\ufffd" in text or re.search(r"\\(?:label|begin\{document|cite)\b", text):
            failures.append(f"{name}: unresolved conversion notation")
        delimiters = re.findall(r"^\$\$$", text, re.M)
        if len(delimiters) % 2 or len(delimiters) != text.count("$$"):
            failures.append(f"{name}: display math delimiters must be paired on separate lines")
        anchors = re.findall(r'<a id="([^"]+)"></a>', text)
        if len(anchors) != len(set(anchors)):
            failures.append(f"{name}: duplicate explicit anchors")
        for target in LINK.findall(text):
            if target.startswith("#") and target[1:] not in anchors:
                failures.append(f"{name}: missing explicit manuscript anchor {target}")
    documents = [ROOT / "README.md", ROOT / "papers/README.md", ROOT / "docs/manuscript-editing.md"]
    documents += manuscripts + sorted((ROOT / "papers").glob("*/README.md"))
    links = 0
    for path in documents:
        for target in LINK.findall(path.read_text(encoding="utf-8")):
            parsed = urlsplit(target)
            if parsed.scheme or not parsed.path:
                continue
            links += 1
            local = (path.parent / unquote(parsed.path)).resolve()
            if not local.is_relative_to(ROOT) or not local.exists():
                failures.append(f"{path.relative_to(ROOT)}: missing local target {target}")
            if local.suffix.lower() == ".pdf":
                failures.append(f"{path.relative_to(ROOT)}: local PDF link {target}")
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    # A staged deletion is absent from the index; an unstaged deletion is absent on disk.
    pdfs = [p for p in tracked if p.lower().endswith(".pdf") and (ROOT / p).exists()]
    if pdfs:
        failures.append(f"Tracked PDF artifacts: {pdfs}")
    for pattern in ("*/content.py", "*/main.tex"):
        for path in (ROOT / "papers").glob(pattern):
            if path.relative_to(ROOT).as_posix() in tracked:
                failures.append(f"Duplicate manuscript source: {path.relative_to(ROOT)}")
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"Passed: {len(manuscripts)} Markdown manuscripts, {links} local links, no tracked PDFs or duplicate prose sources.")


if __name__ == "__main__":
    main()

"""Check the editable manuscript format and local publication links."""
from __future__ import annotations

import json
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r"\]\(([^\s)]+)\)")


def check_math_blocks(text: str, name: Path) -> list[str]:
    """Keep TeX out of the ordinary Markdown escape/emphasis parser."""
    failures = []
    fence = None
    math = False
    for number, line in enumerate(text.splitlines(), 1):
        if (fence is None or math) and r"\operatorname" in line:
            failures.append(f"{name}:{number}: GitHub rejects operatorname; use an upright math name")
        match = re.match(r"^\s*(`{3,}|~{3,})(.*)$", line)
        if match:
            marker, language = match.groups()
            if fence is None:
                fence = marker
                math = language.strip() == "math"
            elif marker[0] == fence[0] and len(marker) >= len(fence) and not language.strip():
                fence = None
                math = False
        elif fence is None and line.lstrip().startswith("$$"):
            failures.append(f"{name}:{number}: use a fenced math block to preserve TeX escapes")
    if fence and math:
        failures.append(f"{name}: unclosed math fence")
    return failures


def main() -> None:
    failures = []
    manuscripts = sorted((ROOT / "docs").glob("*/manuscript.md"))
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    old_paths = [p for p in tracked if p.startswith(("papers/", "provenance/")) or p == "PROVENANCE.md"]
    if old_paths:
        failures.append("Stage the documentation and records moves; old root directories remain tracked.")
    expected = {v["paper"] for v in json.loads((ROOT / "records/markdown-migration.json").read_text(encoding="utf-8"))["manuscripts"]}
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
        anchors = re.findall(r'<a id="([^"]+)"></a>', text)
        if len(anchors) != len(set(anchors)):
            failures.append(f"{name}: duplicate explicit anchors")
        for target in LINK.findall(text):
            if target.startswith("#") and target[1:] not in anchors:
                failures.append(f"{name}: missing explicit manuscript anchor {target}")
    for path in [ROOT / "README.md", *(ROOT / "docs").rglob("*.md")]:
        failures.extend(check_math_blocks(path.read_text(encoding="utf-8"), path.relative_to(ROOT)))
    documents = [ROOT / "README.md", ROOT / "RECORDS.md", ROOT / "docs/README.md", ROOT / "docs/manuscript-editing.md"]
    documents += manuscripts + sorted((ROOT / "docs").glob("*/README.md"))
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
    # A staged deletion is absent from the index; an unstaged deletion is absent on disk.
    pdfs = [p for p in tracked if p.lower().endswith(".pdf") and (ROOT / p).exists()]
    if pdfs:
        failures.append(f"Tracked PDF artifacts: {pdfs}")
    for pattern in ("*/content.py", "*/main.tex"):
        for path in (ROOT / "docs").glob(pattern):
            if path.relative_to(ROOT).as_posix() in tracked:
                failures.append(f"Duplicate manuscript source: {path.relative_to(ROOT)}")
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"Passed: {len(manuscripts)} Markdown manuscripts, {links} local links, no tracked PDFs or duplicate prose sources.")


if __name__ == "__main__":
    main()

"""Explicitly restore original publication artifacts for historical checks."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "provenance/publication-archive.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, help="Exact path or directory prefix in the registry")
    parser.add_argument("--output", type=Path, default=Path("D:/Codex/chaos-publication-archive"))
    parser.add_argument("--materialize", action="store_true", help="Also restore ignored files at their original repository paths")
    args = parser.parse_args()
    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
    prefix = args.prefix.replace("\\", "/").rstrip("/")
    items = [v for v in registry["files"] if v["path"] == prefix or v["path"].startswith(prefix + "/")]
    if not items:
        parser.error("No archived paths match this prefix")
    output = args.output.resolve()
    if output == ROOT or ROOT in output.parents:
        parser.error("The archive output directory must be outside this repository")
    for item in items:
        name = item["path"]
        target = (output / name).resolve()
        if not target.is_relative_to(output):
            raise ValueError(f"Unsafe archive path: {name}")
        data = target.read_bytes() if target.exists() else None
        if data is None or hashlib.sha256(data).hexdigest() != item["sha256"]:
            result = subprocess.run(["git", "cat-file", "blob", item["git_blob"]], cwd=ROOT, capture_output=True)
            if result.returncode == 0:
                data = result.stdout
            else:
                url = f"https://raw.githubusercontent.com/Griboedow/chaos-attractor-lab/{item.get('source_commit', registry['source_commit'])}/{name}"
                with urllib.request.urlopen(url, timeout=60) as response:
                    data = response.read()
            if len(data) != item["bytes"] or hashlib.sha256(data).hexdigest() != item["sha256"]:
                raise ValueError(f"Archive hash mismatch: {name}")
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        if args.materialize:
            local = (ROOT / name).resolve()
            if not local.is_relative_to(ROOT):
                raise ValueError(f"Unsafe repository path: {name}")
            if local.exists() and hashlib.sha256(local.read_bytes()).hexdigest() != item["sha256"]:
                raise FileExistsError(f"Refusing to replace a different local file: {name}")
            local.parent.mkdir(parents=True, exist_ok=True)
            local.write_bytes(data)
    print(f"Verified {len(items)} archived artifacts in {output}")
    if args.materialize:
        print("Restored ignored historical files for legacy checks; do not add them to Git.")


if __name__ == "__main__":
    main()

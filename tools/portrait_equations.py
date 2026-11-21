"""Synchronize the defining equations beside all static portraits in the main README."""
import argparse
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
START = "<!-- portrait-equation "
END = "<!-- /portrait-equation -->"
BLOCK = re.compile(r"\n\n<!-- portrait-equation ([a-z0-9-]+) -->\n.*?\n<!-- /portrait-equation -->", re.S)
PORTRAIT = re.compile(
    r"^!\[[^\]\n]*\]\(docs/assets/portraits/([a-z0-9-]+)\.png\)"
    r"(?:\n\n<!-- drawing-date -->\n[^\n]+\n<!-- /drawing-date -->)?", re.M
)


def render(ident, row):
    source = row["source"]
    links = [("Model and parameters", source)]
    for extra in row.get("sources", []):
        path = Path(extra.split("#", 1)[0])
        case = re.fullmatch(r"case-(\d+)", path.stem)
        if case:
            label = "Case " + case[1]
        elif path.name == "README.md":
            label = path.parent.name.replace("-", " ").capitalize() + " report"
        else:
            label = path.stem.replace("-", " ").replace("_", " ").capitalize()
        links.append((label, extra))
    for _, target in links:
        if not (ROOT / target.split("#", 1)[0]).is_file():
            raise ValueError(f"{ident}: missing equation source {target}")
    if not row.get("text", "").strip():
        raise ValueError(f"{ident}: missing model description")
    parts = [START + ident + " -->"]
    if row.get("tex"):
        tex = row["tex"].strip()
        if "```" in tex:
            raise ValueError(f"{ident}: store the equation without a Markdown fence")
        parts.append("```math\n" + tex + "\n```")
    references = " · ".join(f"[{label}]({target})" for label, target in links)
    parts.append(row["text"].strip() + " " + references)
    return "\n\n".join(parts) + "\n\n" + END


def synchronize(check=False):
    catalog = json.loads((ROOT / "records/portrait-catalog.json").read_text(encoding="utf-8"))["rows"]
    equations = json.loads((ROOT / "records/portrait-equations.json").read_text(encoding="utf-8"))["equations"]
    ids = {r["id"] for r in catalog}
    if set(equations) != ids:
        raise ValueError(f"Equation registry differs from catalog: {sorted(set(equations) ^ ids)}")
    path = ROOT / "README.md"
    source = path.read_text(encoding="utf-8")
    bare = BLOCK.sub("", source)
    if START in bare or END in bare:
        raise ValueError("Malformed equation block in README")
    matches = PORTRAIT.findall(bare)
    if len(matches) != len(ids) or set(matches) != ids:
        raise ValueError("Each cataloged portrait must occur once in README")
    result = PORTRAIT.sub(lambda m: m.group(0) + "\n\n" + render(m[1], equations[m[1]]), bare)
    if check and result != source:
        raise ValueError("Run python tools/portrait_equations.py sync")
    if not check and result != source:
        with path.open("w", encoding="utf-8", newline="\n") as output:
            output.write(result)
    print(f"{len(ids)} portrait equation entries verified")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["sync", "check"])
    synchronize(parser.parse_args().command == "check")

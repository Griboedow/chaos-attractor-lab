"""Group the screened tubes by a full polynomial before selecting representatives."""
import argparse
import json
from functools import lru_cache
from pathlib import Path

from alexander import normalized_braid_polynomial
from search import braid, write_json


@lru_cache(maxsize=8192)
def polynomial(strands, word):
    return tuple(normalized_braid_polynomial(strands,word))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("search_dir", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    screen = json.loads((args.search_dir/"screen.json").read_text(encoding="utf-8"))
    candidates = json.loads((args.search_dir/"survivors.json").read_text(encoding="utf-8"))
    leaders = {}
    for j,case in enumerate(candidates):
        word = braid(case,512)
        if word is None:
            word = braid(case,4096)
        assert word is not None
        poly = polynomial(case["degree"], tuple(word))
        key = case["degree"],poly
        if key not in leaders or leaders[key]["branch_gap_screen"] < case["branch_gap_screen"]:
            leaders[key] = case
        if (j+1) % 5000 == 0:
            print(f"{j+1} tubes, {len(leaders)} degree / polynomial pairs",flush=True)
    screen.update(leaders=list(leaders.values()),classes=len(leaders))
    write_json(args.out,screen)

"""Read the certified node boxes and enumerate the six-sign, one-turn saddle words."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from flint import arb, ctx, fmpq

ROOT = Path(__file__).resolve().parents[4]
PACKET = ROOT / "results/cubic-forced-orbits"
FRAMES = [0, 11, 21, 32, 43, 53]


def multiply(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)] for i in range(2)
    ]


def measure(folder):
    proof = json.loads((folder / "orbits.json").read_text())
    node_file = folder / "shooting-nodes.npz"
    if hashlib.sha256(node_file.read_bytes()).hexdigest() != proof["node_file_sha256"]:
        raise ValueError("Node proposals differ from the orbit certificate")
    ctx.prec = proof["arb_bits"]
    if proof["box_half_width"] != "2^-40" or proof["phase_nodes_per_turn"] != 64:
        raise ValueError("Unexpected node box dimensions")
    radius = arb(2) ** -40
    margin = arb(fmpq(9312, 10000))
    records = []
    with np.load(node_file) as data:
        nodes = data["nodes"]
    for orbit in proof["orbits"]:
        if orbit["type"] not in ("positive saddle", "negative saddle"):
            continue
        word = []
        for k in FRAMES:
            center = nodes[orbit["orbit"], k].real
            x = arb(fmpq(*float(center).as_integer_ratio()), radius)
            if x > margin:
                word.append("+")
            elif x < -margin:
                word.append("-")
            else:
                raise ArithmeticError("Certified sign or separation unresolved")
        changes = sum(a != b for a, b in zip(word, word[1:] + word[:1]))
        sign = (-1) ** (1 + changes // 2)
        if (sign > 0) != (orbit["type"] == "positive saddle"):
            raise ArithmeticError(
                "Expanding multiplier sign disagrees with the measured itinerary"
            )
        records.append(
            {
                "orbit": orbit["orbit"],
                "word": "".join(word),
                "sign_changes": changes,
                "expanding_multiplier_sign": sign,
            }
        )

    # Matrix rows are starting signs in (-,+) order; columns are the ending signs.
    down = [[1, 0], [1, 1]]
    up = [[1, 1], [0, 1]]
    steps = [down, up] * 3
    allowed = set()
    for bits in itertools.product(range(2), repeat=6):
        if all(m[bits[i]][bits[(i + 1) % 6]] for i, m in enumerate(steps)):
            allowed.add("".join("-+"[b] for b in bits))
    product = [[1, 0], [0, 1]]
    for matrix in steps:
        product = multiply(product, matrix)
    if (
        len(allowed) != product[0][0] + product[1][1]
        or len(records) != len(allowed)
        or allowed != {r["word"] for r in records}
    ):
        raise ArithmeticError(
            "The certified saddle words and the six-frame rule differ"
        )
    return {
        "phase_indices": FRAMES,
        "phase_denominator": 64,
        "signs": "sign of Re(z)",
        "certified_abs_real_lower": "0.9312",
        "node_file_sha256": proof["node_file_sha256"],
        "sign_order": ["-", "+"],
        "alternating_allowed_steps": [down, up],
        "one_turn_count_matrix": product,
        "trace": product[0][0] + product[1][1],
        "binary_words_examined": 64,
        "saddle_words": records,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", type=Path, default=PACKET)
    args = parser.parse_args()
    result = measure(args.folder)
    (args.folder / "itineraries.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print("Certified one-turn saddle words:", len(result["saddle_words"]))

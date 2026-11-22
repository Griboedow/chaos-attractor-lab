"""Replay the two seven-switch rows with fixed binary64 inputs and extra precision."""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from flint import arb, ctx


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/art_20261122/memory-transplant"
DONORS = json.loads((RESULTS / "donors.json").read_text(encoding="utf-8"))
RUN = json.loads((RESULTS / "outcomes.json").read_text(encoding="utf-8"))
FIELDS = ("L01", "U04")
PRECISIONS = (160, 320, 640, 1280, 2560)
FIRST, LAST = 62, 82


def replay(job):
    bits, field_id, initial_memory = job
    ctx.prec = bits
    p = {name: arb(value) for name, value in RUN["parameters"].items()}
    row = next(row for row in DONORS["fields"] if row["id"] == field_id)
    x, y = map(arb, row["state"][:2])
    memory = arb(initial_memory)
    assert all(value.is_exact() for value in (*p.values(), x, y, memory))
    slow = 1 - p["relaxation"]
    gain = p["relaxation"] * p["memory_gain"]
    labels = []
    previous = 0
    for first, last in RUN["windows"]:
        low, high = None, None
        for step in range(previous + 1, last + 1):
            intensity = x*x + y*y
            response = intensity / (1 + intensity)
            phase = p["base_detuning"] + p["fast_shift"] * response + memory
            sine, cosine = phase.sin_cos()
            # Take point midpoints every update. These are numerical replays.
            x, y = ((p["pump"] + p["field_retention"] * (x*cosine - y*sine)).mid(),
                    (p["field_retention"] * (x*sine + y*cosine)).mid())
            memory = (slow * memory + gain * response).mid()
            if step >= first:
                if low is None or memory < low:
                    low = memory
                if high is None or memory > high:
                    high = memory
        labels.append("L" if high < RUN["threshold"] else "U" if low > RUN["threshold"] else "?")
        previous = last
    return labels[0] if labels[0] == labels[1] else "?"


def switches(labels):
    return sum(a != b for a, b in zip(labels, labels[1:]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    indices = range(FIRST, LAST + 1)
    keys = [(bits, name, i) for bits in PRECISIONS for name in FIELDS for i in indices]
    # Round the dose on the binary64 grid once, as in the original experiment.
    jobs = [(bits, name, 5.5 + 3.7*i/192) for bits, name, i in keys]
    donor_fields = [row["id"] for row in DONORS["fields"]]
    swaps = [(160, name, donor["state"][2]) for name in donor_fields for donor in DONORS["fields"]]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        output = dict(zip(keys, pool.map(replay, jobs)))
        swap_labels = "".join(pool.map(replay, swaps))
    records = []
    for name in FIELDS:
        row = next(row for row in RUN["rows"] if row["field"] == name)
        labels = {str(b): "".join(output[b, name, i] for i in indices) for b in PRECISIONS}
        records.append({
            "field": name,
            "reference_binary64": row["doses"][FIRST:LAST + 1],
            "numpy_binary64": row["numpy_doses"][FIRST:LAST + 1],
            "labels_by_precision_bits": labels,
            "switches_by_precision_bits": {b: switches(s) for b, s in labels.items() if "?" not in s},
        })
    summary = {
        "case": RUN["case"], "inclusive_zero_based_indices": [FIRST, LAST],
        "dose_binary64_formula": "5.5 + (3.7 * j) / 192",
        "dose_binary64_hex": [(5.5 + 3.7*i/192).hex() for i in indices],
        "windows": RUN["windows"], "threshold": RUN["threshold"], "rows": records,
        "donor_swap_precision_bits": 160, "donor_order": donor_fields,
        "donor_swap_rows": [swap_labels[i:i+len(donor_fields)]
                            for i in range(0, len(swap_labels), len(donor_fields))],
    }
    (RESULTS / "precision.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

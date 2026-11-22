"""Replay the donor interventions using NumPy complex arithmetic."""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
donors = json.loads((RESULTS / "art_20261122/memory-transplant/donors.json").read_text())
settings = json.loads((RESULTS / "wow_search_20261122/orbit/case-86780.json").read_text())
p = settings["parameters"]
states = np.array([entry["state"] for entry in donors["fields"]])
fields = states[:, 0] + 1j * states[:, 1]
doses = [5.5 + 3.7 * j / 192 for j in range(193)]


def classify(starting_memories):
    light = np.repeat(fields, len(starting_memories))
    memory = np.tile(starting_memories, len(fields))
    step = 0
    labels = []
    for first, last in ((6145, 8192), (24577, 32768)):
        minimum = np.full(light.size, np.inf)
        maximum = np.full(light.size, -np.inf)
        while step < last:
            step += 1
            intensity = light.real**2 + light.imag**2
            response = intensity / (1 + intensity)
            phase = p["base_detuning"] + p["fast_shift"] * response + memory
            light = p["pump"] + p["field_retention"] * light * np.exp(1j * phase)
            memory = (1 - p["relaxation"]) * memory + p["relaxation"] * p["memory_gain"] * response
            if step >= first:
                minimum = np.minimum(minimum, memory)
                maximum = np.maximum(maximum, memory)
        labels.append(np.where(maximum < 7.5, "L", np.where(minimum > 7.5, "U", "?")))
    same = labels[0] == labels[1]
    result = np.where(same, labels[1], "?").reshape(len(fields), len(starting_memories))
    return ["".join(row) for row in result]


if __name__ == "__main__":
    print(json.dumps({
        "numpy": np.__version__,
        "doses": classify(doses),
        "swaps": classify(states[:, 2])
    }))

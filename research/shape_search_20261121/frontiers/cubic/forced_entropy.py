"""Read the three-sink certificate and compute exact fixed-point and entropy bounds."""

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def multiply(left, right):
    return [
        [sum(left[i][k] * right[k][j] for k in range(2)) for j in range(2)]
        for i in range(2)
    ]


def braid_matrix(word):
    generators = {
        1: [[1, 1], [0, 1]],
        -1: [[1, -1], [0, 1]],
        2: [[1, 0], [-1, 1]],
        -2: [[1, 0], [1, 1]],
    }
    result = [[1, 0], [0, 1]]
    for letter in word:
        result = multiply(generators[letter], result)
    return result


def bounds(packet):
    low, high = packet["records"]
    assert packet["case"] == 48000001
    assert low["gain"] == [3, 100] and high["gain"] == [17, 50]
    for record in (low, high):
        assert record["checks"][-1]["success"]
    assert low["dyadic_crossing_word"] == []
    word = high["dyadic_crossing_word"]
    assert word == [-1, 2] * 3
    matrix = braid_matrix(word)
    assert matrix == [[5, -8], [-8, 13]]
    trace = matrix[0][0] + matrix[1][1]
    last, current = 2, trace
    counts = []
    for turns in range(1, 7):
        counts.append(
            {
                "turns": turns,
                "trace": current,
                "finite_fixed_points_at_least": current + 3,
            }
        )
        last, current = current, trace * current - last
    d = Fraction(1, 32)
    lower = 9 * (Fraction(1, 2) - d) / 5
    assert lower > 6 * 3 * d
    assert 1200 * d * lower > 30 + Fraction(1, 4)
    return {
        "case": packet["case"],
        "gain": high["gain"],
        "word": word,
        "matrix": matrix,
        "entropy_lower_exact": "log(9 + 4 sqrt(5))",
        "entropy_lower_approx": math.acosh(trace / 2),
        "fixed_points_of_iterates": counts,
        "family_sufficient_gain_per_harmonic_index": 1200,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("finite_gains", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = bounds(json.loads(args.finite_gains.read_text(encoding="utf-8")))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

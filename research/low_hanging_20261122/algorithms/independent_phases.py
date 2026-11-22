"""Recheck the six no-return phases using the archived normal-vector projectors."""

from itertools import product
import json
from pathlib import Path
import sys

import mpmath as mp

ROOT = Path(__file__).resolve().parents[3]
ARCHIVE = ROOT / "results/raar-symbolic-dynamics"
sys.path.insert(0, str(ARCHIVE / "source"))
from certify_raar_coverings import Model, exactfloat, matmul, matvec


def run():
    mp.iv.dps = 60
    iv = mp.iv
    data = json.loads((ARCHIVE / "runs/raar-covering-refined01/input.json").read_bytes())
    model = Model(data["params"])
    hull = [iv.mpf([a, b]) for a, b in [("-0.02", "0.08"), ("-2.61", "-2.38"), ("4.02", "4.19")]]
    gap = iv.mpf("0.01")
    phases = 0
    for c, b, r in zip(data["centers"], data["bases"], data["radii"]):
        center = [exactfloat(x, iv) for x in c]
        matrix = [[exactfloat(row[j], iv) * exactfloat(r[j], iv) for j in range(3)] for row in b]
        wide = matvec(matrix, [iv.mpf([-1, 1])] * 3)
        assert all((center[k] + wide[k]) in hull[k] for k in range(3))
        widths = [4, 8, 2]
        offsets = [iv.mpf([-1 / n, 1 / n]) for n in widths]
        for index in product(*(range(n) for n in widths)):
            middle = [iv.mpf(-1 + (2 * i + 1) / n) for i, n in zip(index, widths)]
            point = [x + y for x, y in zip(center, matvec(matrix, middle))]
            derivative = matrix
            for _ in range(6):
                box = [x + y for x, y in zip(point, matvec(derivative, offsets))]
                _, jacobian, _ = model.step(box)
                point, _, _ = model.step(point, False)
                derivative = matmul(jacobian, derivative)
                image = [x + y for x, y in zip(point, matvec(derivative, offsets))]
                assert any(y.a > (z + gap).b or y.b < (z - gap).a for y, z in zip(image, hull))
                phases += 1
    print(f"Independent 60-decimal check passed for {phases} cell/phase combinations.")


if __name__ == "__main__":
    run()

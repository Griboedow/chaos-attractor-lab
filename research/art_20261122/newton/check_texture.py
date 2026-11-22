"""Compare the torus classifier to direct Arb capture at full, unfolded scales."""
import json
import random

from sectors import OUT, Q, acb, arb, capture, load_inputs, poles
from torus import inputs, texture


def main():
    rng = random.Random(6902915)
    a, roots = load_inputs()
    af, rf = inputs()
    records = {"seed": 6902915, "point_count": 0, "certified_matches": 0,
               "unresolved": 0, "wrong": 0, "matches_by_label": [0]*15}
    for idx in range(256):
        if idx % 2:
            angle = -3.141592653589793-.01 + rng.uniform(-.0425, .0425)
        else:
            angle = rng.uniform(-3.141592653589793, 3.141592653589793)
        frac = rng.random()
        for cell in range(9):
            scale = frac+cell
            label, _ = texture(angle, scale, af, rf)
            # No period folding here. For every |1/t| >= 1000, the external
            # recurrence bounds |Psi(t)-1/t| by 6e-12, uniformly in scale.
            w = 1000*Q()**(-arb(scale))*acb(0, -arb(angle)).exp()
            outcome = capture(w, arb("6.1e-12"), a, roots)
            records["point_count"] += 1
            if outcome is None:
                records["unresolved"] += 1
                continue
            ref = outcome[0]
            if ref >= 6:
                ref = 6+(ref-6-1)%9
            if ref == label:
                records["certified_matches"] += 1
                records["matches_by_label"][ref] += 1
            else:
                records["wrong"] += 1
    assert records["wrong"] == 0
    assert min(records["matches_by_label"]) > 0
    sectors = json.loads((OUT/"sectors.json").read_text())["sectors"]
    native_checks = 0
    for sector in sectors:
        lo, hi = sector["direction_degrees_closed"]
        for angle in [lo, (lo+hi)//2, hi]:
            direction = acb(0, angle*arb.pi()/180).exp()
            for radius in ["1e-14", "1e-15", "1e-16"]:
                offset = arb(radius)*direction
                for pole, expected in zip(poles(a), sector["root_labels_1based"]):
                    result = capture(pole+offset, arb("1e-44"), a, roots)
                    assert result is not None and result[0]+1 == expected
                    native_checks += 1
    records["direct_common_displacement_checks"] = native_checks
    (OUT/"texture-check.json").write_text(json.dumps(records, indent=2)+"\n")
    print(records)


if __name__ == "__main__":
    main()

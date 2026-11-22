"""Certify the five outcomes obtained by translating all Newton poles together."""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from math import comb
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "research/attractor_followup_2/newton"))
from radial import Q, acb, arb, capture, load_inputs

OUT = ROOT / "results/art_20261122/newton"
WINDOWS = [(149, 153), (207, 211)]


def poles(a):
    principal = (a/6)**(arb(1)/5)
    return [principal * acb(0, 2*arb.pi()*j/5).exp() for j in range(5)]


def regular_bound(pole, residue):
    # Desingularize the pole before enclosing the binomial remainder.
    h = arb("0.001")
    p = abs(pole)
    d0 = abs(30*pole**4)
    dvar = sum(6*comb(5, j)*p**(5-j)*h**(j-1) for j in range(2, 6))
    nvar = sum(5*comb(6, j)*p**(6-j)*h**(j-1) for j in range(1, 7))
    assert d0 > dvar
    return (nvar+abs(residue)*dvar/h)/(d0-dvar)


def shell_bounds(ps, a):
    for p in ps:
        residue = p*(5*p-6)/30
        assert regular_bound(p, residue) < 2
        assert abs(residue)*Q()**120/1000 < arb("0.001")
        assert abs(residue)*Q()**120/1000 > arb("1.79568e-14")
    floor = arb(989)
    assert abs(a)*(Q()+1/floor)/((6-abs(a)/floor**5)*floor**4) < arb("1e-12")
    assert 10*Q()**120 + 6*arb("1e-12") < arb("1e-8")


def cover(task):
    pole_id, window, cell, cells = task
    a, roots = load_inputs()
    p = poles(a)[pole_id]
    residue = p*(5*p-6)/30
    lo_angle, hi_angle = WINDOWS[window]
    half_angle = (hi_angle-lo_angle)*arb.pi()/(360*cells)
    angle = lo_angle*arb.pi()/180 + (2*cell+1)*half_angle
    unit = residue/abs(residue)*acb(0, -angle).exp()
    stack = [(arb(1000), arb(1200), 0)]
    label = None
    leaves = 0
    max_depth = 0
    while stack:
        lo, hi, depth = stack.pop()
        mid = (lo+hi)/2
        radius = (hi-lo)/2 + hi*half_angle + arb("1e-8")
        outcome = capture(unit*mid, radius, a, roots)
        if outcome is None:
            if depth >= 18:
                raise ValueError(f"Unclassified tile at {pole_id=} {window=} {cell=} {mid=}")
            stack.extend([(lo, mid, depth+1), (mid, hi, depth+1)])
            continue
        assert outcome[0] < 6 and label in (None, outcome[0]), (task, label, outcome)
        label = outcome[0]
        max_depth = max(max_depth, depth)
        leaves += 1
    return pole_id, window, label, leaves, max_depth


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    a, roots = load_inputs()
    ps = poles(a)
    shell_bounds(ps, a)
    cells = 192
    jobs = [(p, w, c, cells) for p in range(5) for w in range(2) for c in range(cells)]
    if args.workers > 1:
        with ProcessPoolExecutor(args.workers) as pool:
            parts = list(pool.map(cover, jobs, chunksize=8))
    else:
        parts = list(map(cover, jobs))
    rows = []
    for w, window in enumerate(WINDOWS):
        labels = []
        counts = []
        for p in range(5):
            rows_p = [r for r in parts if r[0:2] == (p, w)]
            label, = {r[2] for r in rows_p}
            labels.append(label+1)
            counts.append(sum(r[3] for r in rows_p))
        assert len(set(labels)) == 5
        rows.append({"direction_degrees_closed": window, "root_labels_1based": labels,
                     "tiles_at_each_pole": counts})
    result = {"case": 6902, "common_radius": "1.79568e-14", "n_min": 120,
              "reference_R": [1000, 1200], "all_poles_regular_part_bound": 2,
              "reference_landing_error_bound": "1e-8", "angle_cells_per_sector": cells,
              "sectors": rows}
    dest = OUT / "sectors.json"
    if args.verify:
        assert json.loads(dest.read_text()) == json.loads(json.dumps(result))
    else:
        dest.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(rows))


if __name__ == "__main__":
    main()

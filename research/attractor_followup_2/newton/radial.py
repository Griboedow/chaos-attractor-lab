"""Certify every basin on each of infinitely many radial intervals at a Newton pole."""
import argparse
import json
from math import comb
from pathlib import Path

from flint import acb, arb, ctx

ROOT = Path(__file__).resolve().parents[3]
PRIOR = ROOT / "results/wow_search_20261122/boundary/certificate.json"
CASE_DIR = ROOT / "results/attractor_followup_2/newton"
HORIZONTAL = {
    "east": [10472, 2797, 1954, 7911, 5899, 1864, 1851, 8433, 8485, 1911,
             7028, 4876, 6402, 3233, 8003],
    "west": [6042, 2321, 10200, 2001, 209, 2024, 227, 11474, 1236, 1285,
             798, 500, 11374, 11404, 11057],
}
AXES = {"east": 1, "west": -1, "north": 1j, "south": -1j}


def Q():
    return arb(5) / 6


def center(record):
    a, b = record["center_dyadic"]
    return acb(arb((int(a[0]), a[1])), arb((int(b[0]), b[1])))


def square(z, radius):
    r = arb(0, radius.upper())
    return z + acb(r, r)


def value(z, a):
    if abs(z) > 3:
        return Q()*z + a*(Q()*z-1)/(6*z**5-a)
    return (5*z**6-a)/(6*z**5-a)


def derivative(z, a):
    if abs(z) > 3:
        d = 6*z**5-a
        return Q() + a*(Q()*d - (Q()*z-1)*30*z**4)/d**2
    return 30*z**4*(z**6-a*z+a)/(6*z**5-a)**2


def capture(start, radius, a, roots):
    z = start.mid()
    rad = (radius + abs(z-start)).upper()
    for step in range(250):
        if rad > 110:
            return None
        if abs(z).lower() < 2:
            if abs(z)+rad < arb("1e-7"):
                return 6+(-step)%9, step
            for j, root in enumerate(roots):
                if abs(z-root)+rad < arb("9.99e-8"):
                    return j, step
        sup = abs(derivative(square(z, rad), a)).upper()
        if not sup.is_finite():
            return None
        y = value(z, a)
        z = y.mid()
        rad = (rad*sup + abs(y-z)).upper()
    return None


def load_inputs():
    ctx.prec = 240
    original = json.loads(PRIOR.read_bytes())
    a = square(center(original["parameter_certificate"]), arb("1e-142"))
    roots = [center(rec) for rec in original["root_certificates"]]
    return a, roots



def bounds(a, pole, residue):
    h = arb("0.001")
    p = abs(pole)
    d0 = abs(30*pole**4)
    dvar = sum(6*comb(5, j)*p**(5-j)*h**(j-1) for j in range(2, 6))
    nvar = sum(5*comb(6, j)*p**(6-j)*h**(j-1) for j in range(1, 7))
    assert d0 > dvar
    regular = (nvar + abs(residue)*dvar/h) / (d0-dvar)
    assert regular < 1
    floor = arb(989)
    tail = abs(a)*(Q()+1/floor) / ((6-abs(a)/floor**5)*floor**4)
    assert tail < arb("1e-12")
    assert abs(residue)*Q()**120/1000 < h
    assert 10*Q()**120 + 6*arb("1e-12") < arb("1e-8")
    return {
        "pole_regular_part_upper": 1,
        "large_z_linear_remainder_upper": "1e-12",
        "landing_error_upper": "1e-8",
        "outer_radius_n120_approx": float(abs(residue)*Q()**120/1000),
    }


def horizontal(direction, unit, a, roots):
    ray = unit * acb(AXES[direction]).conjugate()
    rows = []
    # Simultaneously cover a radial interval and a small in-plane angular change.
    assert arb("1e-7") + 1200*arb("5e-10") + arb("1e-8") < arb("1e-6")
    for j in HORIZONTAL[direction]:
        ref = ray*(1000 + arb(j)/60)
        outcome = capture(ref, arb("1e-6"), a, roots)
        assert outcome is not None, (direction, j)
        label, steps = outcome
        rows.append({"j": j, "landing_label": label, "steps_to_capture": steps})
    assert sorted(row["landing_label"] for row in rows) == list(range(15))
    return sorted(rows, key=lambda row: row["landing_label"])


def sector(job):
    direction, cell, cells = job
    a, roots = load_inputs()
    pole = (a/6)**(arb(1)/5)
    residue = pole*(5*pole-6)/30
    # Subdivide +/-12 degrees exactly, using pi in radians.
    half_sector = arb.pi()/15
    half_cell = half_sector/cells
    theta = -half_sector + (2*cell+1)*half_cell
    unit = residue/abs(residue)*acb(AXES[direction]).conjugate()
    ray = unit*acb(0, -theta).exp()
    stack = [(arb(1000), arb(1200), 0)]
    leaves = 0
    max_depth = 0
    max_steps = 0
    label = 1 if direction == "north" else 4
    while stack:
        lo, hi, depth = stack.pop()
        mid = (lo+hi)/2
        rad = (hi-lo)/2 + hi*half_cell + arb("1e-8")
        outcome = capture(ray*mid, rad, a, roots)
        if outcome is None:
            assert depth < 18, (direction, cell, float(mid))
            stack.extend([(lo, mid, depth+1), (mid, hi, depth+1)])
            continue
        assert outcome[0] == label, (direction, cell, outcome)
        leaves += 1
        max_depth = max(max_depth, depth)
        max_steps = max(max_steps, outcome[1])
    return direction, leaves, max_depth, max_steps


def direct_control(a, roots, pole, residue):
    # Finite control in the original plane, prior to any normal-form reduction.
    checked = 0
    for direction, js in HORIZONTAL.items():
        for j in js:
            start_radius = 1000 + arb(j)/60
            native_unit = acb(AXES[direction])
            ref = residue/abs(residue) * native_unit.conjugate() * start_radius
            ref_label, _ = capture(ref, arb("1e-6"), a, roots)
            for n in [120, 121, 250]:
                z = pole + native_unit*abs(residue)*Q()**n/start_radius
                # After n+1 exact-formula steps the label should match the shell.
                for _ in range(n+1):
                    z = (5*z**6-a)/(6*z**5-a)
                    z = z.mid()
                assert abs(z-ref) < arb("1e-8")
                outcome = capture(z, arb("1e-15"), a, roots)
                assert outcome is not None and outcome[0] == ref_label
                checked += 1
    return checked


def main():
    from concurrent.futures import ProcessPoolExecutor

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    a, roots = load_inputs()
    pole = (a/6)**(arb(1)/5)
    residue = pole*(5*pole-6)/30
    results = {
        "case": 6902, "n_min": 120, "R_bounds": [1000, 1200],
        "horizontal_R_halfwidth": "1e-7", "horizontal_angular_halfwidth_rad": "5e-10",
        "vertical_angular_halfwidth_degrees": 12,
        "shell_estimates": bounds(a, pole, residue),
        "horizontal": {key: horizontal(key, residue/abs(residue), a, roots)
                       for key in HORIZONTAL},
        "vertical": {},
        "direct_original_plane_samples": direct_control(a, roots, pole, residue),
    }
    cells = 512
    jobs = [(name, cell, cells) for name in ("north", "south") for cell in range(cells)]
    if args.workers == 1:
        parts = list(map(sector, jobs))
    else:
        with ProcessPoolExecutor(args.workers) as pool:
            parts = list(pool.map(sector, jobs, chunksize=8))
    for name, root in [("north", 1), ("south", 4)]:
        tiles = [row for row in parts if row[0] == name]
        results["vertical"][name] = {
            "root_index_0based": root, "angle_cells": cells,
            "certified_tiles": sum(row[1] for row in tiles),
            "deepest_radial_bisection": max(row[2] for row in tiles),
            "max_steps_to_capture": max(row[3] for row in tiles),
        }
    target = CASE_DIR / "certificate.json"
    if args.verify:
        assert json.loads(target.read_text()) == results
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results["vertical"]))


if __name__ == "__main__":
    main()

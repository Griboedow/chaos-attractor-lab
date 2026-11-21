"""Sample forward iterates in the fixed Euclidean embedding of the solid torus."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numba import njit

from search import write_json


@njit(cache=True)
def trajectory(p, amplitude, offset, modes, scale, rho, twist, phase, initial, burn, steps):
    theta, x, y = initial
    xyz = np.empty((steps, 3))
    growth = 0.
    for j in range(steps+burn):
        cx, cy = 0., 0.
        for k,u,v in modes:
            s,c = math.sin(k*theta),math.cos(k*theta)
            cx += u*c - v*s
            cy += u*s + v*c
        s,c = math.sin(twist*theta+phase),math.cos(twist*theta+phase)
        x,y = cx/scale + rho*(c*x-s*y), cy/scale + rho*(s*x+c*y)
        if j >= burn:
            growth += math.log(p*(1+amplitude*math.cos(p*theta)))
        theta = (p*theta+amplitude*math.sin(p*theta)+offset) % (2*math.pi)
        if j >= burn:
            xyz[j-burn] = ((1.15+x)*math.cos(theta), (1.15+x)*math.sin(theta), y)
    return xyz, growth/steps


def sample(row, initial, steps):
    return trajectory(row["degree"], row["base_ripple"], row["base_offset"],
                      np.array(row["harmonics"], dtype=np.int64),row["scale"],
                      row["fiber_radius"],row["fiber_turns"],row["fiber_phase"],
                      initial,1000,steps)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificates", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--steps", type=int, default=60000)
    ap.add_argument("--cases", nargs="+", type=int)
    args = ap.parse_args()
    cases = json.loads(args.certificates.read_text(encoding="utf-8"))["rows"]
    if args.cases is not None:
        cases = [r for r in cases if r["case"] in args.cases]
        assert len(cases) == len(set(args.cases))
    args.out.mkdir(parents=True, exist_ok=True)
    initial_states = [(0.71,0.24,-0.03),(4.09,-0.12,0.65),(2.23,0.01,0.84)]
    results = []
    for row in cases:
        stats = []
        for k, initial in enumerate(initial_states):
            points, rate = sample(row,initial,args.steps)
            assert np.isfinite(points).all()
            if k == 0:
                path = args.out / f"case-{row['case']}.npz"
                np.savez_compressed(path, points=points)
            stats.append({"initial": list(initial), "base_exponent": rate,
                          "mean": points.mean(axis=0).tolist(),
                          "covariance_spectrum": np.linalg.eigvalsh(np.cov(points.T)).tolist()})
        results.append({"family":"Fourier disk-bundle map","case":row["case"],
                        "degree":row["degree"],"determinant":row["verified_core_determinant"],
                        "cloud":path.name, "runs":stats})
    write_json(args.out/"manifest.json", {"steps_per_initial":args.steps, "burn":1000, "rows":results})
    print(f"{len(cases)} geometries sampled from {len(initial_states)} initial states")

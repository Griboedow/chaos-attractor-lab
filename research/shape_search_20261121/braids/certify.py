"""Certify a Fourier core braid by rational bounds and an exact polygonal diagram."""
import argparse
import json
import math
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path

from search import color_determinant, write_json

Q = 10**15


def atan_bounds(inverse):
    partial = sum(F((-1)**k, (2*k+1) * inverse**(2*k+1)) for k in range(20))
    return partial, partial + F(1, 41 * inverse**41)


@lru_cache(maxsize=None)
def roots(n):
    # Machin's formula, with alternating-series bounds at both arguments.
    a, b = atan_bounds(5), atan_bounds(239)
    lo, hi = 2*(16*a[0]-4*b[1])/n, 2*(16*a[1]-4*b[0])/n
    assert 0 < lo <= hi < F(1,20)
    trig = []
    for parity in (0,1):
        last = 10 + parity
        poly = sum((-1)**j * lo**(2*j+parity) / math.factorial(2*j+parity)
                   for j in range(6))
        error = hi-lo + hi**(last+1)/math.factorial(last+1)
        assert error < F(1, 4*Q)
        trig.append(round(poly*Q))
    # Each rounded component differs by < 1/Q from the true step.
    # Iteration uses a unit exact root. Rounded product error < 2/Q,
    # stored modulus <= 2, so 10/Q per step bounds accumulated error.
    c, s = trig
    x, y = Q, 0
    points = []
    for _ in range(n):
        assert x*x + y*y <= 4*Q*Q
        points.append((x,y))
        x,y = (x*c-y*s)//Q, (x*s+y*c)//Q
    return points


def mesh(case, nodes):
    n = case["degree"]*nodes
    unit = roots(n)
    harmonic = case["harmonics"]
    denom = Q * case["scale"]
    amplitude = sum(abs(x)+abs(y) for k,x,y in harmonic)
    speed = sum(abs(k)*(abs(x)+abs(y)) for k,x,y in harmonic)
    vertices = []
    for j in range(n):
        # A small rational perturbation avoids exactly aligned diagram vertices.
        real, imag = (j*1171 + j*j*17) % 2039, (j*769 + j*j*53) % 2029
        for k, a, b in harmonic:
            x,y = unit[(k*j) % n]
            real += a*x-b*y
            imag += a*y+b*x
        vertices.append((real, imag))
    vertex_error = F(10*n*amplitude+4096, denom)
    # A linear interpolant differs from a Lipschitz function by <= L*h/2.
    error = vertex_error + F(22*speed, 7*case["scale"]*n)
    rad = max(x*x+y*y for x,y in vertices)
    radius_upper = F(math.isqrt(rad)+1, denom) + error
    gap_square = min((x-vertices[(j+nodes*k)%n][0])**2 +
                     (y-vertices[(j+nodes*k)%n][1])**2
                     for j,(x,y) in enumerate(vertices)
                     for k in range(1,case["degree"]//2+1))
    gap_lower = F(math.isqrt(gap_square),denom) - 2*error
    rho = F(str(case["fiber_radius"]))
    assert 0 < rho and radius_upper + rho < 1 and gap_lower > 2*rho
    assert gap_lower > 2*error and radius_upper + error < 1
    return vertices, radius_upper, gap_lower, error


def exact_word(vertices, degree, nodes, shear):
    a,b = shear
    diagram = [(a*x+b*y, a*y-b*x) for x,y in vertices]
    n = len(vertices)
    word = []
    for j in range(nodes):
        left = [diagram[j + k*nodes] for k in range(degree)]
        right = [diagram[(j+1 + k*nodes)%n] for k in range(degree)]
        events = []
        for u,v in combinations(range(degree),2):
            x0, x1 = left[u][0]-left[v][0], right[u][0]-right[v][0]
            if x0*x1 == 0:
                return None
            if x0*x1 < 0:
                phase = F(x0, x0-x1)
                events.append((phase, u, v))
        order = sorted(range(degree), key=lambda i:left[i][0])
        events.sort()
        for idx, (phase,u,v) in enumerate(events):
            if idx and events[idx-1][0] == phase:
                return None
            pos = min(order.index(u),order.index(v))
            assert set(order[pos:pos+2]) == {u,v}
            u,v = order[pos:pos+2]
            h0,h1 = left[u][1]-left[v][1],right[u][1]-right[v][1]
            height = h0*(phase.denominator-phase.numerator) + h1*phase.numerator
            assert height
            word.append(pos+1 if height>0 else -pos-1)
            order[pos:pos+2] = [v,u]
        assert order == sorted(range(degree), key=lambda i:right[i][0])
    return word


def certify(case, nodes):
    vertices, radius, gap, error = mesh(case,nodes)
    for shear in ((101,29),(503,173),(151,67),(2111,157)):
        word = exact_word(vertices,case["degree"],nodes,shear)
        if word is not None:
            break
    else:
        raise ValueError("No generic projection in the configured list")
    determinant = color_determinant(case["degree"], word)
    # Distinct numerical screens may alias a close crossing. The exact diagram wins.
    return {**case, "nodes_per_strand": nodes, "shear": list(shear),
            "radius_upper_ppm": math.ceil(radius*10**6),
            "branch_gap_lower_ppm": math.floor(gap*10**6),
            "polygon_error_upper_ppm": math.ceil(error*10**6),
            "exact_word": word, "verified_core_determinant": determinant}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("screen", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--nodes", type=int, default=512)
    args = parser.parse_args()
    data = json.loads(args.screen.read_text(encoding="utf-8"))
    rows = [certify(case, args.nodes) for case in data["leaders"]]
    rows.sort(key=lambda r: (r["degree"],r["verified_core_determinant"],r["case"]))
    write_json(args.out, {"draws": data["draws"], "seed": data["seed"],
               "screen_passes": data["survivors"], "screen_unresolved": data["unresolved_crossings"],
               "rows": rows})
    determinants = {r["verified_core_determinant"] for r in rows}
    print(f"{len(rows)} embedded trapping tubes, {len(determinants)} distinct first-core determinants")

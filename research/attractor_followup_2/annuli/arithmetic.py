"""Check arithmetic nonrepetition, regularity and disjointness for Case 13745."""

import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TRAPS = ROOT / "results/wow_search_20261122/symmetry/annuli.json"


def order_two_positive(n):
    assert n > 0
    return (n & -n).bit_length() - 1


def gaussian_order(x, y):
    """Order at 1+i in Q(i), from the ordinary 2-order of the Gaussian norm."""
    norm = x*x + y*y
    return order_two_positive(norm.numerator) - order_two_positive(norm.denominator)


def dominant(terms):
    """Term orders are intercept + degree*v. Compare at v=-1 and below."""
    intercept, degree = terms[-1]
    assert all(degree > d and intercept-degree < c-d for c, d in terms[:-1])


def residue(x):
    return x.numerator * pow(x.denominator, -1, 5) % 5


def multiply(z, w):
    a, b = z
    c, d = w
    return ((a*c-b*d) % 5, (a*d+b*c) % 5)


def power(z, n):
    result = (1, 0)
    for _ in range(n):
        result = multiply(result, z)
    return result


def next_mod5(z, terms):
    conj = (z[0], -z[1])
    x, y = 0, 0
    for a, b, i, j in terms:
        dx, dy = multiply((residue(a), residue(b)), multiply(power(z, i), power(conj, j)))
        x, y = (x+dx) % 5, (y+dy) % 5
    return x, y


def check():
    case = json.loads(TRAPS.read_text())["13745"]
    n, lam, alpha, beta, gamma, omega, delta, e, h = [
        Fraction(float.fromhex(x)) for x in case["binary64_hex"]]
    assert n == 17
    terms = [(lam, omega, 1, 0), (alpha, delta, 2, 1),
             (beta/2, 0, 18, 0), (beta/2, 0, 1, 17),
             (gamma, 0, 0, 16), (e, h, 0, 33)]
    orders = [(gaussian_order(x, y), a+b) for x, y, a, b in terms]
    dominant(orders)
    assert orders[-1] == (-110, 33)
    # Compare every term of both Wirtinger derivatives at once. The last
    # conjugate derivative term is the sole minimum and dominates both sums.
    dz = [(gaussian_order(x*a, y*a), a+b-1) for x, y, a, b in terms if a]
    db = [(gaussian_order(x*b, y*b), a+b-1) for x, y, a, b in terms if b]
    dominant(dz + db)
    assert db[-1] == (-110, 32)
    starts = [tuple(map(Fraction, entry["start"])) for entry in case["entrance_tests"]]
    vals = [gaussian_order(*z) for z in starts]
    assert max(vals) <= -1 and 33*max(vals)-110 < min(vals)
    print(f"Initial (1+i)-orders: {vals}. Recurrence v -> 33*v - 110, nonsingular derivative.")
    # Distinct levels of time have disjoint ranges of orders. At equal times
    # the affine recurrence separates starts with different initial orders.
    for i, j in combinations(range(len(starts)), 2):
        if vals[i] != vals[j]:
            continue
        pair = tuple(tuple(map(residue, starts[k])) for k in (i, j))
        visited = set()
        while pair not in visited:
            assert pair[0] != pair[1]
            visited.add(pair)
            pair = tuple(next_mod5(z, terms) for z in pair)
        print(f"Starts {i+1} and {j+1}: distinct at all equal times ({len(visited)} mod-5 pairs).")


if __name__ == "__main__":
    check()

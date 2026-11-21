"""Compute normalized Alexander polynomials and exhaust torus-knot polynomials of each span."""
import argparse
import json
import math
from pathlib import Path

from flint import fmpz_poly as Poly

from search import write_json

T = Poly([0,1])
ONE = Poly([1])


def determinant(matrix):
    matrix = [row.copy() for row in matrix]
    sign,previous = 1,ONE
    for k in range(len(matrix)-1):
        pivot = next((i for i in range(k,len(matrix)) if matrix[i][k]), None)
        if pivot is None:
            return Poly([])
        if pivot != k:
            matrix[k],matrix[pivot] = matrix[pivot],matrix[k]
            sign = -sign
        for i in range(k+1,len(matrix)):
            for j in range(k+1,len(matrix)):
                cross = matrix[k][k]*matrix[i][j] - matrix[i][k]*matrix[k][j]
                matrix[i][j],remainder = divmod(cross,previous)
                assert not remainder
        previous = matrix[k][k]
    return sign*matrix[-1][-1]


def normalized_braid_polynomial(p,word):
    matrix = [[Poly([int(i==j)]) for j in range(p)] for i in range(p)]
    unit = ONE
    for g in word:
        i = abs(g)-1
        l,r = matrix[i],matrix[i+1]
        if g > 0:
            matrix[i] = [(1-T)*a+T*b for a,b in zip(l,r)]
            matrix[i+1] = l
        else:
            # Keep polynomial entries by clearing one Laurent denominator.
            matrix = [[T*v for v in row] for row in matrix]
            unit = T*unit
            matrix[i] = [T*b for b in r]
            matrix[i+1] = [a+(T-1)*b for a,b in zip(l,r)]
    for i in range(p):
        matrix[i][i] -= unit
    poly = determinant([row[:-1] for row in matrix[:-1]])
    coeff = [int(c) for c in poly.coeffs()]
    first = next(j for j,c in enumerate(coeff) if c)
    coeff = coeff[first:]
    if coeff[0] < 0:
        coeff = [-v for v in coeff]
    assert coeff == coeff[::-1] and abs(sum(coeff)) == 1
    return coeff


def torus_matches(coeff):
    span = len(coeff)-1
    if span == 0:
        return [[1,1]]
    matches = []
    for a in range(2,math.isqrt(span)+2):
        if span % (a-1):
            continue
        b = 1 + span//(a-1)
        if b < a or math.gcd(a,b) != 1:
            continue
        candidate,rem = divmod((T**(a*b)-1)*(T-1),(T**a-1)*(T**b-1))
        assert not rem
        if candidate == Poly(coeff):
            matches.append([a,b])
    return matches


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificates", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    assert normalized_braid_polynomial(2,[1]*3) == [1,-1,1]
    assert normalized_braid_polynomial(3,[1,-2]*2) == [1,-3,1]
    assert normalized_braid_polynomial(3,[1,1,1,-2]) == [1,-1,1]
    assert torus_matches([1,-1,1]) == [[2,3]] and not torus_matches([1,-3,1])
    data = json.loads(args.certificates.read_text(encoding="utf-8"))
    rows = []
    for case in data["rows"]:
        poly = normalized_braid_polynomial(case["degree"],case["exact_word"])
        det = abs(sum(v*(-1)**j for j,v in enumerate(poly)))
        assert det == case["verified_core_determinant"]
        rows.append({"case":case["case"],"degree":case["degree"],"determinant":det,
                     "alexander_low_to_high":poly,"torus_polynomial_matches":torus_matches(poly)})
    distinct = {tuple(r["alexander_low_to_high"]) for r in rows}
    excluded = {tuple(r["alexander_low_to_high"]) for r in rows if not r["torus_polynomial_matches"]}
    write_json(args.out,{"distinct_polynomials":len(distinct),
               "distinct_excluding_torus_knots":len(excluded), "rows":rows})
    print(f"{len(distinct)} distinct core polynomials, {len(excluded)} exclude all torus knots")

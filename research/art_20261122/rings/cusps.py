"""Certify a Whitney cusp of the 17-fold polynomial map and its capture inside K."""
import json
import sys
from pathlib import Path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "research/low_hanging_20261122/annuli"))
import verify as base

PACKET = ROOT / "results/art_20261122/rings"


def differentiate(terms, direction):
    derivative = {}
    for (a, b), c in terms.items():
        if a:
            key = (a-1, b)
            derivative[key] = derivative.get(key, 0) + a*direction*c
        if b:
            key = (a, b-1)
            derivative[key] = derivative.get(key, 0) + b*direction.conjugate()*c
    return derivative


def cross(v, w):
    return (v.conjugate()*w).imag


def cusp_system(terms, imag):
    polynomials = {}
    for a, b in ((1, 0), (0, 1), (2, 0), (1, 1), (0, 2),
                 (3, 0), (2, 1), (1, 2), (0, 3)):
        if a:
            polynomials[a, b] = differentiate(polynomials.get((a-1, b), terms), 1)
        else:
            polynomials[a, b] = differentiate(polynomials.get((a, b-1), terms), imag)

    def evaluate(z, screening=False):
        values = (sum(c*z**a*z.conjugate()**b for (a, b), c in p.items())
                  for p in polynomials.values())
        dx, dy, xx, xy, yy = [next(values) for _ in range(5)]
        J = cross(dx, dy)
        Jx = cross(xx, dy) + cross(dx, xy)
        Jy = cross(xy, dy) + cross(dx, yy)
        if screening:
            return J, dx*Jy - dy*Jx
        xxx, xxy, xyy, yyy = values
        Jxx = cross(xxx, dy) + 2*cross(xx, xy) + cross(dx, xxy)
        Jxy = cross(xxy, dy) + cross(xx, yy) + cross(dx, xyy)
        Jyy = cross(xyy, dy) + 2*cross(xy, yy) + cross(dx, yyy)
        H = dx.real*Jy - dy.real*Jx
        Hx = xx.real*Jy + dx.real*Jxy - xy.real*Jx - dy.real*Jxx
        Hy = xy.real*Jy + dx.real*Jyy - yy.real*Jx - dy.real*Jxy
        return J, H, [[Jx, Jy], [Hx, Hy]], dx.real
    return evaluate


def inputs():
    case = json.loads(base.TRAPS.read_text())["13745"]
    coeff = [base.exact(float.fromhex(x)) for x in case["binary64_hex"][1:]]
    l, a, b, g, w, d, e, h = coeff
    terms = {(1, 0): acb(l, w), (2, 1): acb(a, d), (18, 0): b/2,
             (1, 17): b/2, (0, 16): g, (0, 33): acb(e, h)}
    return coeff, terms


def count_cusps(evaluate, center, known, jacobian):
    """Exclude extra singularities in one closed sector of both trapping annuli."""
    edges = json.loads(base.TRAPS.read_text())["13745"]["annuli"]
    slope = (arb.pi()/17).tan()
    total = 0
    for lo, hi in edges:
        lower, upper = base.exact(lo), base.exact(hi)
        x0 = lower*(arb.pi()/17).cos()
        y1 = upper*(arb.pi()/17).sin()
        boxes = [(x0.union(upper), (-y1).union(y1))]
        while boxes:
            x, y = boxes.pop()
            z = acb(x, y)
            total += 1
            if abs(z) < lower or abs(z) > upper or abs(y) > slope*x:
                continue
            if abs(x-center.real) < known and abs(y-center.imag) < known:
                continue
            if abs(jacobian(z).det()) > 0:
                continue
            J, tangent = evaluate(z, screening=True)
            if abs(J) > 0 or abs(tangent.real) > 0 or abs(tangent.imag) > 0:
                continue
            assert max(x.rad().upper(), y.rad().upper()) > base.exact("1e-7")
            if x.rad() > y.rad():
                boxes.extend(((x.lower().union(x.mid()), y), (x.mid().union(x.upper()), y)))
            else:
                boxes.extend(((x, y.lower().union(y.mid())), (x, y.mid().union(y.upper()))))
    return total


def verify(build, bits):
    ctx.prec = bits
    coeff, terms = inputs()
    evaluate = cusp_system(terms, base.IMAG)
    value, jacobian = build(coeff)
    w = json.loads((PACKET / "cusp-witness.json").read_text())
    center = acb(*[base.exact(float.fromhex(x)) for x in w["source_binary64"]])
    cusp_width = base.exact("1e-15")
    box = base.square(center, cusp_width)

    def equations(c, r):
        J, H, _, _ = evaluate(c)
        return acb(J, H), arb_mat(evaluate(base.square(c, r))[2])

    base.inverse_branch(equations, center, cusp_width, acb(0), arb(0), 1)
    known = base.exact("1e-4")
    base.inverse_branch(equations, center, known, acb(0), arb(0), 1)
    cells = count_cusps(evaluate, center, known, jacobian)
    _, _, derivatives, Px = evaluate(box)
    assert Px < -2 and arb_mat(derivatives).det() > 95000
    target, width = center, base.exact(w["target_half_width"])
    assert width > cusp_width
    ratios = []
    for steps, coords, extent, tiles in w["backward_covers"]:
        source = acb(*[base.exact(float.fromhex(x)) for x in coords])
        radius = base.exact(extent)

        def image(c, r):
            return base.iterate(c, r, steps, value, jacobian)

        ratio, _ = base.inverse_branch(image, source, radius, target, width, tiles)
        ratios.append(ratio)
        target, width = source, radius
    initial = json.loads((base.PACKET / "witness.json").read_text())
    rotation = (base.IMAG * 2*arb.pi()*w["rotation"]/17).exp()
    offset = target/rotation - acb(*map(base.exact, initial["centers"][0]))
    assert max(abs(offset.real).upper(), abs(offset.imag).upper()) + arb(2).sqrt()*width < (
        base.exact(initial["half_width_in_billionths"][0])/10**9)
    cusp_image = value(box)
    assert abs(cusp_image) > base.exact("0.41127") and abs(cusp_image) < base.exact("0.41129")
    print(f"{build.__name__}, {bits} bits: ordinary cusp at {complex(center)}, "
          f"image {complex(cusp_image)}, {cells} cells exclude additional cusps in the sector; "
          f"largest inverse inclusion ratio {float(max(ratios)):.6f}.")


if __name__ == "__main__":
    verify(base.polynomial, 192)
    verify(base.expanded_polynomial, 256)

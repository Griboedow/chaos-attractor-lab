"""Check two inverse branches in the trapped outer ring with python-flint 0.9.0."""

import json
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, arb_mat, ctx


ROOT = Path(__file__).resolve().parents[3]
PACKET = ROOT / "results/low_hanging_20261122/annuli"
TRAPS = ROOT / "results/wow_search_20261122/symmetry/annuli.json"
IDENTITY = arb_mat([[1, 0], [0, 1]])
IMAG = acb(0, 1)


def exact(value):
    rational = Fraction(value)
    return arb(rational.numerator) / rational.denominator


def square(center, radius):
    delta = arb(0, radius.upper())
    return center + acb(delta, delta)


def row_norm(matrix, i):
    return abs(matrix[i, 0]) + abs(matrix[i, 1])


def norm(matrix):
    return max(row_norm(matrix, i).upper() for i in range(2))


def jac_columns(dx, dy):
    return arb_mat([[dx.real, dy.real], [dx.imag, dy.imag]])


def polynomial(coefficients):
    lam, alpha, beta, gamma, omega, delta, e, h = coefficients
    L, C, E = acb(lam, omega), acb(alpha, delta), acb(e, h)

    def value(z):
        u = z.real**2 + z.imag**2
        return (L + C*u + beta*(z**17).real)*z + gamma*z.conjugate()**16 + E*z.conjugate()**33

    def jacobian(z):
        u = z.real**2 + z.imag**2
        z16, z17 = z**16, z**17
        A = L + 2*C*u + beta*z17.real + (17*beta/2)*z17
        B = C*z*z + (17*beta/2)*z*z16.conjugate() + 16*gamma*z.conjugate()**15 + 33*E*z.conjugate()**32
        return jac_columns(A+B, IMAG*(A-B))

    return value, jacobian


def expanded_polynomial(coefficients):
    """Independent value and real-direction derivatives from the expanded monomials."""
    lam, alpha, beta, gamma, omega, delta, e, h = coefficients
    terms = [(acb(lam, omega), 1, 0), (acb(alpha, delta), 2, 1),
             (beta/2, 18, 0), (beta/2, 1, 17), (gamma, 0, 16), (acb(e, h), 0, 33)]

    def value(z):
        bar = z.conjugate()
        return sum(c*z**a*bar**b for c, a, b in terms)

    def jacobian(z):
        bar = z.conjugate()
        dz = sum(c*a*z**(a-1)*bar**b for c, a, b in terms if a)
        db = sum(c*b*z**a*bar**(b-1) for c, a, b in terms if b)
        return jac_columns(dz+db, IMAG*(dz-db))

    return value, jacobian


def iterate(center, radius, steps, value, jacobian):
    # The integral mean-value formula bounds images of the original square.
    # The running derivative encloses D(F^j) throughout that square.
    ref, derivative = center, IDENTITY
    for _ in range(steps):
        dx = arb(0, (radius*row_norm(derivative, 0)).upper())
        dy = arb(0, (radius*row_norm(derivative, 1)).upper())
        derivative = jacobian(ref + acb(dx, dy)) * derivative
        ref = value(ref)
    return ref, derivative


def inverse_branch(H, source, sr, target, tr, tiles):
    # Verify a preimage in source + [-sr,sr]^2 for every point of the target square.
    image, linear = H(source, arb(0))
    inverse = linear.inv()
    A = arb_mat([[inverse[i, j].mid() for j in range(2)] for i in range(2)])
    assert abs(A.det()) > 0
    residual = image - target
    displacement = A * arb_mat([[residual.real], [residual.imag]])
    offset = max(abs(displacement[i, 0]).upper() for i in range(2))
    tile = sr/tiles
    q = arb(0)
    for i in range(tiles):
        for j in range(tiles):
            center = source + acb(tile*(2*i+1-tiles), tile*(2*j+1-tiles))
            _, derivative = H(center, tile)
            q = max(q, norm(IDENTITY-A*derivative))
    assert q < 1
    inclusion = (offset + norm(A)*tr)/sr + q
    assert inclusion < 1
    return inclusion.upper(), (norm(A)/(1-q)).upper()


def verify(build, bits):
    ctx.prec = bits
    spec = json.loads((PACKET / "witness.json").read_text())
    trap = json.loads(TRAPS.read_text())[str(spec["case"])]
    n, *coef = [exact(float.fromhex(x)) for x in trap["binary64_hex"]]
    assert n == 17 and spec["return_steps"] == 18 and spec["rotation_sectors"] == 1
    value, jacobian = build(coef)
    rotation = (IMAG * (2*arb.pi()/17)).exp()
    rotation_matrix = jac_columns(rotation, IMAG*rotation)

    def H(center, radius):
        image, derivative = iterate(center, radius, 18, value, jacobian)
        return rotation*image, rotation_matrix*derivative

    centers = [acb(exact(x), exact(y)) for x, y in spec["centers"]]
    radii = [exact(x)/10**9 for x in spec["half_width_in_billionths"]]
    assert len(centers) == len(radii) == len(spec["grid_per_axis"]) == 10
    lower, upper = map(exact, trap["annuli"][1])
    first = centers[0]
    total = arb(1)
    ratios = []
    for i, (center, radius, tiles) in enumerate(zip(centers, radii, spec["grid_per_axis"])):
        position = abs(square(center, radius))
        assert position > lower and position < upper
        previous = max(i-1, 0)
        ratio, lipschitz = inverse_branch(H, center, radius, centers[previous], radii[previous], tiles)
        ratios.append(ratio)
        if i == 0:
            assert lipschitz < exact("0.60")
        else:
            total *= lipschitz

    assert max(ratios) < exact("0.90") and total < exact("0.035")
    assert abs((centers[1]-first).real) > radii[1]+radii[0]
    displacement = centers[-1]-first
    assert max(abs(displacement.real).upper(), abs(displacement.imag).upper()) + radii[-1] < radii[0]

    # Neighboring centers are the closest pair of rotations. Enclose each
    # rotated square inside its circumscribed Euclidean disk.
    nearest = 2*abs(first)*(arb.pi()/17).sin()
    assert nearest > exact("0.28") and nearest > 2*arb(2).sqrt()*radii[0]

    # Localize the fixed point of the primary inverse branch. Its first and
    # twelfth F-images rule out the remaining proper divisors of 306.
    small = exact("1e-12")
    _, inverse_lip = inverse_branch(H, first, small, first, small, 1)
    assert inverse_lip < 1
    _, M = H(first, small)
    trace = M[0, 0]+M[1, 1]
    disc = (trace**2-4*M.det()).sqrt()
    positive, negative = (trace+disc)/2, (trace-disc)/2
    assert positive > exact("4.2980") and positive < exact("4.2982")
    assert negative > exact("-3.6754") and negative < exact("-3.6752")
    assert abs(jacobian(square(first, small)).det()) > 0
    twelfth, J = iterate(first, small, 12, value, jacobian)
    twelfth_domain = square(twelfth, small*norm(J))
    assert abs(twelfth_domain) < exact("0.74") and abs(square(first, small)) > exact("0.77")
    print(f"{build.__name__}, {bits} bits: all 10 inverse covers verified.")
    print(f"Largest inclusion ratio {float(max(ratios)):.6f}, alternate branch Lipschitz bound {float(total):.6f}.")
    print("Least period 306; both components have interior; entropy >= log(2)/162.")


if __name__ == "__main__":
    verify(polynomial, 160)
    verify(expanded_polynomial, 256)

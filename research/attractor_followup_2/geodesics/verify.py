"""Bound transverse expansion over a closed straight geodesic of winding (7, 7, 1)."""
import argparse
import json
from math import ceil, floor, log10, pi, sqrt
from pathlib import Path

import numpy as np
from flint import arb, ctx
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/attractor_followup_2/geodesics"
POWERS = (2, 7, 12, 24, 36)
CELLS = 16384


def ellipse(power):
    c = arb(3)/8
    b = 1 + c - arb(2)**(-power)
    p = (c*c - 1 - b*b)/2
    h = ((b-p)*(b+p)).sqrt()
    r = -1-p
    trace = (1+r)**2 + (p-r)**2 + (7*r)**2 + 54*h*h
    area = arb(99).sqrt()*h
    ratio_sum = trace/area
    return (ratio_sum + (ratio_sum**2-4).sqrt())/2


def advance(matrix, rate, length):
    # Exact constant-q propagator for q=rate^2 in the (J,J') basis.
    a, b, c, d = matrix
    angle = rate*length
    ch, sh = angle.cosh(), angle.sinh()
    return (ch*a+sh*c/rate, ch*b+sh*d/rate,
            sh*rate*a+ch*c, sh*rate*b+ch*d)


def enclose(power):
    ctx.prec = 144
    k = ellipse(power)
    quarter = arb.pi()/2
    low = high = (arb(1), arb(0), arb(0), arb(1))
    left = arb(0)
    rate_left = k
    for i in range(1, CELLS+1):
        right = quarter*(arb(i)/CELLS)**4
        rate_right = k/(1+(k*k-1)*right.sin()**2)
        dx = right-left
        # The coefficient decreases on the entire quarter period.
        low = advance(low, rate_right.lower(), dx)
        high = advance(high, rate_left.upper(), dx)
        left, rate_left = right, rate_right
    def full_log10(matrix):
        half_trace = 1+2*matrix[1]*matrix[2]
        return 14*half_trace.acosh()/arb(10).log()
    lower, upper = full_log10(low), full_log10(high)
    # Round a second time to printable bounds and assert strict inequalities
    # with Arb itself, avoiding any dependence on binary64 formatting.
    lo = floor(float(lower.lower())*100)-1
    hi = ceil(float(upper.upper())*100)+1
    assert lower > arb(lo)/100 and upper < arb(hi)/100
    return lo, hi, float(k.mid())


def independent(power):
    """Integrate through a full half-period using the surface's implicit Hessian."""
    c = .375
    b = 1+c-2.0**(-power)
    p = (c*c-1-b*b)/2
    # Heron form avoids cancellation near the degenerate triangle.
    h2 = (1+b+c)*(1+b-c)*(1+c-b)*(b+c-1)/4
    h = sqrt(h2)
    r = -1-p
    U = np.array([1+r, p-r, 7*r])
    V = np.array([-h, 2*h, -7*h])
    gram = np.array([[U@U, U@V], [U@V, V@V]])
    eigenvalues, axes = np.linalg.eigh(gram)
    peak = np.arctan2(axes[0, 0], axes[1, 0])
    width = sqrt(eigenvalues[0]/eigenvalues[1])
    direction = np.array([1., 1., 1./7])
    wave = np.array([1., -1., 7.])
    def jacobi(angle, state):
        x = angle+peak
        sn, cs = np.sin(x), np.cos(x)
        # Evaluate the gradient and Hessian-vector product of the three waves
        # directly. Each wave phase increases at unit speed along direction.
        grad = -(np.array([sn, p*sn+h*cs, 0.]) + (r*sn-h*cs)*wave)
        hess_direction = -(np.array([cs, p*cs-h*sn, 0.]) + (r*cs+h*sn)*wave)
        cross = np.cross(direction, grad)
        transverse = cross / np.linalg.norm(cross)
        # The normal curvature along the straight ray is zero, so the
        # mixed second fundamental form alone gives -K times squared speed.
        q = (transverse@hess_direction)**2 / (grad@grad)
        a, b0, c0, d0 = state
        return [c0, d0, q*a, q*b0]
    state = np.array([1., 0., 0., 1.])
    # Force the independent solver to sample both narrow ends of this interval.
    points = np.unique(np.clip([0, 8*width, pi/2, pi-8*width, pi], 0, pi))
    for start, end in zip(points, points[1:]):
        sol = solve_ivp(jacobi, (start, end), state, method="DOP853",
                        rtol=1e-10, atol=1e-11, max_step=min(.08, (end-start)/16))
        if not sol.success:
            raise ArithmeticError(sol.message)
        state = sol.y[:, -1]
    return float(14*np.arccosh((state[0]+state[3])/2)/np.log(10))


def plot(rows, dest):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a,b) = plt.subplots(1, 2, figsize=(10.5, 3.8), layout="constrained")
    powers = [row["distance_power"] for row in rows]
    values = [row["cartesian_log10"] for row in rows]
    a.plot(powers, values, "o-", c="#19777c", lw=2.1, label="Jacobi equation")
    a.plot(powers, [14*log10(pi*pi*row["ellipse_ratio"]) for row in rows],
           "--", c="#bd672e", label="Sharp limiting law")
    a.set(xlabel="m, with second weight 11/8 − 2⁻ᵐ", ylabel="log₁₀ of the one-lap multiplier",
          title="Same lap length, growing instability")
    a.legend(frameon=False)
    x=np.linspace(-.3,.3,1601)
    for row,color in zip(rows[:3],["#36496d","#bd672e","#19777c"]):
        k=row["ellipse_ratio"]
        rate=k/(1+(k*k-1)*np.sin(x)**2)
        b.semilogy(x, rate**2, color=color, label=f"m = {row['distance_power']}")
    b.set(xlabel="Phase offset from the curvature peak (radians)", ylabel="Jacobi coefficient",
          title="Where the expansion builds up")
    b.legend(frameon=False)
    for ax in (a,b):
        ax.spines[["top","right"]].set_visible(False)
    fig.savefig(dest,dpi=155)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Write the compact enclosures")
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()
    rows = []
    for power in POWERS:
        low, high, ratio = enclose(power)
        control = independent(power)
        assert low/100 < control < high/100, (power, low, control, high)
        rows.append(dict(distance_power=power, lower_hundredths=low, upper_hundredths=high,
                         ellipse_ratio=ratio, cartesian_log10=control))
        print(f"2^-{power}: {low/100:.2f} < log10(multiplier) < {high/100:.2f}; independent {control:.4f}")
    assert rows[-1]["lower_hundredths"] > 9000
    path = OUT / "bounds.json"
    if args.write:
        path.write_text(json.dumps(rows, indent=2)+"\n", encoding="utf8", newline="\n")
    else:
        for actual, saved in zip(rows, json.loads(path.read_text()), strict=True):
            for field in ("distance_power","lower_hundredths","upper_hundredths"):
                assert actual[field] == saved[field]
    if args.figure:
        plot(rows,args.figure)


if __name__ == "__main__":
    main()

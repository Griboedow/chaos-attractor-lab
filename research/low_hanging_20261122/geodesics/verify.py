"""Count straight closed geodesics and independently integrate one bent return."""
import argparse
import json
from collections import Counter
from fractions import Fraction
from math import gcd, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_bvp, solve_ivp

N = 7
C = 3 / 8
B = 9 / 8


def cycles(step, sheets):
    """Enumerate the deck-translation orbits, independently of the gcd count."""
    unused = set(range(sheets))
    wraps = []
    while unused:
        sheet = min(unused)
        count = 0
        while sheet in unused:
            unused.remove(sheet)
            count += 1
            sheet = (sheet + step) % sheets
        wraps.append(count)
    assert len(wraps) == gcd(sheets, step)
    return wraps


def length_counts(n, equal_weights):
    """Integer squared lengths in units of (2 pi)^2."""
    result = Counter()
    for slope_y in (-1, 1):
        for slope_w in ((-1, 0, 1) if equal_weights else (-1, 1)):
            step = slope_w - 1 + slope_y
            for turns in cycles(step, n):
                dz = Fraction(step * turns, n)
                length_squared = 2 * turns**2 + dz**2
                assert length_squared.denominator == 1
                result[int(length_squared)] += 2  # two reflected phase choices
    return dict(sorted(result.items()))


def chart(x, h, hp):
    """Chart with y = x + pi + h and w in (-pi, 0)."""
    y = x + pi + h
    u = 1 + hp
    cos_w = -(np.cos(x) + B * np.cos(y)) / C
    if np.any(abs(cos_w) >= 1):
        raise ValueError("Periodic solution left the chosen smooth chart")
    w = -np.arccos(cos_w)
    sin_w = np.sin(w)
    gx = -np.sin(x) - C * sin_w
    gy = -B * np.sin(y) + C * sin_w
    gz = -N * C * sin_w
    zp = -(gx + gy * u) / gz
    z = (w - x + y) / N
    quadratic = -np.cos(x) - B * np.cos(y) * u**2 - C * cos_w * (1 - u + N * zp)**2
    normal = -quadratic / (gx**2 + gy**2 + gz**2)
    return y, z, zp, normal * (gy - gx * u)


def periodic_state(nodes, tolerance):
    def equation(x, state):
        h, hp = state
        return np.vstack((hp, chart(x, h, hp)[3]))
    mesh = np.linspace(0, 2*pi, nodes)
    fit = solve_bvp(equation, lambda left, right: left-right, mesh,
                    np.zeros((2, len(mesh))), tol=tolerance, max_nodes=4096)
    if not fit.success:
        raise RuntimeError(fit.message)
    return fit


def cartesian(t, state):
    """Independent affine-time, unit-speed geodesic equation in the ambient space."""
    x, y, z, vx, vy, vz = state
    w = x - y + N*z
    gradient = np.array([-np.sin(x) - C*np.sin(w),
                         -B*np.sin(y) + C*np.sin(w), -N*C*np.sin(w)])
    hessian_vv = -np.cos(x)*vx**2 - B*np.cos(y)*vy**2 - C*np.cos(w)*(vx-vy+N*vz)**2
    return np.r_[state[3:], -hessian_vv * gradient / (gradient @ gradient)]


def measure():
    x = np.linspace(0, 2*pi, 4097)
    fit = periodic_state(97, 1e-10)
    h, hp = fit.sol(x)
    y, z, zp, _ = chart(x, h, hp)
    length = simpson(np.sqrt(1 + (1+hp)**2 + zp**2), x=x)
    initial_velocity = np.array([1, 1+hp[0], zp[0]])
    initial = np.r_[[0, y[0], z[0]], initial_velocity / np.linalg.norm(initial_velocity)]
    run = solve_ivp(cartesian, (0, length), initial, method="DOP853",
                    atol=2e-13, rtol=2e-13, max_step=0.10)
    if not run.success:
        raise RuntimeError(run.message)
    error = np.linalg.norm(run.y[:, -1] - initial - np.array([2*pi, 2*pi, 0, 0, 0, 0]))
    repeat = periodic_state(65, 1e-8).sol(x)
    assert np.max(abs(repeat[0]-h)) < 5e-9 and error < 1e-8
    lower = 2*pi*sqrt(2)
    upper = 2*pi*sqrt(2 + (1/(3*N))**2)
    assert lower < length < upper
    assert length_counts(7, True) == {2: 14, 99: 6, 102: 2, 107: 2}
    assert length_counts(7, False) == {99: 6, 107: 2}
    for n in (5, 7, 11, 13, 17, 19, 23, 25):
        assert sum(length_counts(n, True).values()) == 2*n+10
        assert sum(length_counts(n, False).values()) == 8
        assert min(length_counts(n, False)) == 2*n*n+1
    receipt = {
        "n": N, "b": "9/8", "c": "3/8",
        "straight_length_squared_in_full_turn_units": {
            "b=1": length_counts(N, True), "b=9/8": length_counts(N, False)},
        "shortest_straight_length": 2*pi*sqrt(99),
        "curved_class_length_bounds": [lower, upper],
        "measured_geodesic_length": length,
        "initial_xyz_unit_velocity": initial.tolist(),
        "return_error_position_and_velocity": error,
        "coarse_fine_max_difference_y": float(np.max(abs(repeat[0]-h))),
        "measured_height_excursion": float(np.ptp(z)),
    }
    return receipt, x, z


def draw(path, x, z):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
    for weights, counts, color, offset in (
        ("Second weight = 1", length_counts(N, True), "#19777c", -0.17),
        ("Second weight = 9/8", length_counts(N, False), "#cb6a2a", 0.17),
    ):
        keys = (2, 99, 102, 107)
        bars = a.bar(np.arange(4)+offset, [counts.get(k, 0) for k in keys],
                     0.32, color=color, label=weights)
        a.bar_label(bars, labels=[str(counts[k]) if k in counts else "" for k in keys], fontsize=9)
    a.set(xticks=range(4), xticklabels=[f"{2*pi*sqrt(k):.2f}" for k in keys],
          xlabel="Primitive lap length", ylabel="Straight closed orbits",
          title="Sixteen straight routes disappear", ylim=(0,17))
    a.legend(frameon=False, fontsize=8)
    b.axhline(0, color="#19777c", lw=2, label="Original straight orbit")
    b.plot(x/(2*pi), z-pi/(2*N), color="#cb6a2a", lw=2.6, label="Bent geodesic, numerical")
    b.set(xlabel="Turns in x", ylabel="Height relative to the original orbit",
          title="A short return survives", xlim=(0,1))
    b.legend(frameon=False, fontsize=8)
    for ax in (a,b):
        ax.spines[["top","right"]].set_visible(False)
    fig.savefig(path, dpi=155)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    result, x, z = measure()
    report = json.dumps(result, indent=2) + "\n"
    if args.receipt:
        args.receipt.write_text(report, encoding="utf-8", newline="\n")
    if args.figure:
        draw(args.figure, x, z)
    print(report)

"""Follow the period-three saddle that collides with the stationary cubic site."""

import argparse
import json
from pathlib import Path

import numpy as np
from flint import arb
from scipy.integrate import solve_ivp
from scipy.optimize import root

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "results/cubic-elliptic-sites"
T = 2 * np.pi


def integrate(gain, initial, tight=False):
    def field(t, h):
        forcing = 0.25 * np.sin(t) - 4 * np.sin(3 * t) + 3j * np.sin(6 * t)
        return -gain * (6 * h + 3 * h**2).conj() * (forcing - 2 + 3 * h**2 + h**3)

    precision = 2.4e-14 if tight else 7e-13
    seeds = np.atleast_1d(initial)
    solution = solve_ivp(
        field,
        (0, 3 * T),
        seeds,
        t_eval=T * np.arange(1, 4),
        method="DOP853",
        rtol=precision,
        atol=abs(seeds) * precision * 0.01,
        max_step=0.07 if tight else 0.16,
    )
    if not solution.success:
        raise ArithmeticError(solution.message)
    return solution.y[0] if np.ndim(initial) == 0 else solution.y


def linear_trace(gain):
    def rhs(t, y):
        a = gain * (12 - 1.5 * np.sin(t) + 24 * np.sin(3 * t))
        b = -18 * gain * np.sin(6 * t)
        return (np.array([[a, b], [b, -a]]) @ y.reshape(2, 2)).ravel()

    sol = solve_ivp(
        rhs, (0, T), np.eye(2).ravel(), method="DOP853", rtol=3e-14, atol=2e-15
    )
    if not sol.success:
        raise ArithmeticError(sol.message)
    return np.trace(sol.y[:, -1].reshape(2, 2))


def experiment(packet):
    row = packet["resonances"][1]
    center = float(row["gain_center"])
    C = complex(*(float(arb(s)) for s in row["coefficient"]))
    slope, height = (float(arb(s)) for s in row["frame_second_column"])

    def physical(v):
        return v.real + (slope + 1j * height) * v.imag

    def internal(h):
        return h.real - slope * h.imag / height + 1j * h.imag / height

    angle = (np.angle(C) + np.pi / 2) / 3
    records = []
    for dx in (-5e-7, -2e-7, -8e-8, -3.2e-8, 3.2e-8, 8e-8, 2e-7, 5e-7):
        gain = center + dx
        delta = 3 * np.arccos(linear_trace(gain) / 2) - 2 * np.pi
        scale = abs(delta) / abs(C)
        guess = np.sign(delta) * np.exp(1j * angle)

        def mismatch(pair, gain=gain, scale=scale, delta=delta):
            z = complex(*pair)
            f = internal(integrate(gain, physical(z * scale))[-1]) / scale
            e = (f - z) / abs(delta)
            return [e.real, e.imag]

        fit = root(mismatch, [guess.real, guess.imag], tol=2e-7, options={"eps": 2e-8})
        z = complex(*fit.x) * scale
        orbit = internal(integrate(gain, physical(z), tight=True))
        relative = abs(orbit[-1] - z) / scale
        if abs(z) < 0.4 * scale or relative > 1e-5:
            raise ArithmeticError(f"Unresolved orbit at {gain}: {relative}")
        phases = np.array([z, *orbit[:2]])
        separation = min(abs(phases - np.roll(phases, 1))) / scale
        if separation < 1.2:
            raise ArithmeticError("Lost the three distinct phases")
        multipliers = []
        for step in (2e-4 * scale, 1e-4 * scale):
            probes = physical(z + step * np.array([1, -1, 1j, -1j]))
            answers = internal(integrate(gain, probes, tight=True)[:, -1])
            columns = (answers[::2] - answers[1::2]) / (2 * step)
            values = np.linalg.eigvals(np.array([columns.real, columns.imag]))
            if not np.all(np.isreal(values)) or not min(values) < 1 < max(values):
                raise ArithmeticError(
                    "The three-period orbit was not numerically a saddle"
                )
            multipliers.append(sorted(float(x) for x in values))
        records.append(
            {
                "gain_offset": dx,
                "gain": gain,
                "linear_detuning": delta,
                "asymptotic_radius": scale,
                "three_turn_multipliers": multipliers,
                "normalized_radii": (abs(phases) / scale).tolist(),
                "relative_closure": relative,
                "scaled_minimum_separation": separation,
                "phase_zero_points": [[v.real, v.imag] for v in phases],
                "physical_offsets": [[v.real, v.imag] for v in physical(phases)],
            }
        )
        print(dx, abs(phases) / scale, relative, flush=True)
    return {
        "gain_center": row["gain_center"],
        "resonant_coefficient": [C.real, C.imag],
        "frame_second_column": [slope, height],
        "records": records,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", type=Path, default=OUT / "resonances.json")
    parser.add_argument("--out", type=Path, default=OUT / "three-cycle.json")
    args = parser.parse_args()
    data = experiment(json.loads(args.packet.read_text()))
    args.out.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

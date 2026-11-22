"""Recover periodic itineraries of the motor rows and check their full-state return."""

import argparse
import json
from itertools import permutations
from pathlib import Path

import numpy as np
from check_motor_asymmetry import augmented
from motor_rows import source_force_constraint_field
from motor_symmetry_passages import (
    LARGE,
    LEVELS,
    TRAIN,
    distances,
    field,
    sample,
    section,
    templates,
    visits,
)
from scipy.integrate import solve_ivp


def vector(t, x, p):
    out = np.empty_like(x)
    field(x, p, out)
    return out


def section_event(t, x, p):
    return x[:4].mean()


section_event.direction = 1


def get_cycle(states, p):
    times, points = section(states, p, 0.0, 0.1)
    # Demand agreement on 100 consecutive crossings in the full-state metric.
    lag = next(
        l
        for l in range(1, 161)
        if np.max(np.linalg.norm(points[-100:] - points[-100 - l : -l], axis=1))
        < 0.0002
    )
    period = times[-1] - times[-1 - lag]
    sol = solve_ivp(
        vector,
        (0, period * 4.5),
        states[-1],
        args=(p,),
        events=section_event,
        method="DOP853",
        rtol=2e-12,
        atol=2e-14,
        max_step=0.04,
    )
    assert sol.success
    t, x = sol.t_events[0], sol.y_events[0]
    start = x[-1 - lag].copy()
    period = t[-1] - t[-1 - lag]

    # One Newton correction for the phase-fixed shooting problem, repeated until small.
    for iteration in range(8):
        z = np.r_[start, np.eye(11).ravel(), 0.0]
        sol = solve_ivp(
            lambda t, y: augmented(y, p, 4),
            (0, period),
            z,
            method="DOP853",
            rtol=8e-13,
            atol=8e-15,
            max_step=0.04,
        )
        assert sol.success
        end = sol.y[:11, -1]
        monodromy = sol.y[11:-1, -1].reshape(11, 11)
        residual = np.r_[end - start, start[:4].mean()]
        if np.max(np.abs(residual)) < 2e-12:
            break
        derivative = np.zeros((12, 12))
        derivative[:11, :11] = monodromy - np.eye(11)
        derivative[:11, -1] = vector(0, end, p)
        derivative[-1, :4] = 0.25
        correction = np.linalg.solve(derivative, -residual)
        start += correction[:11]
        period += correction[-1]
    assert np.max(np.abs(residual)) < 1e-9
    eig = np.abs(np.linalg.eigvals(monodromy))
    neutral = np.argmin(np.abs(eig - 1))
    return (
        start,
        float(period),
        {
            "section_lag": lag,
            "shooting_residual_inf": float(np.max(np.abs(residual))),
            "flow_multiplier_modulus": float(eig[neutral]),
            "leading_transverse_multiplier_modulus": float(
                np.max(np.delete(eig, neutral))
            ),
        },
    )


def labels(d, ratio):
    ranks = np.argsort(d, axis=1)
    near, second = np.take_along_axis(d, ranks[:, :2], axis=1).T
    valid = (near < 0.075) & (second > ratio * near)
    seq = ranks[valid, 0]
    if not len(seq):
        return seq
    return seq[np.r_[True, np.diff(seq) != 0]]


def act(x, perm):
    full = np.r_[x, -x[8:].sum()]
    ix = np.array(perm)
    return full[np.r_[ix, ix + 4, ix + 8]][:11]


def screen(out, p, source, independent, trees):
    records = []
    for reverse in (False, True):
        rates = np.linspace(1.5, 1.54, 17)
        if reverse:
            rates = rates[::-1]
        start = independent.copy() if reverse else source[-1].copy()
        for b in rates:
            p[2] = b
            states = sample(start, p, 400.0, 5000.0, 0.01)
            start = states[-1]
            time, points = section(states, p, 0.0, 0.1)
            result = visits(distances(trees, points), time, 0.075, 2.0)
            records.append(
                {
                    "direction": "down" if reverse else "up",
                    "b": b,
                    "visits": result["visits"],
                    "changes": result["changes"],
                    "edge": result["edge"],
                    "face": result["face"],
                    "opposite": result["opposite"],
                }
            )
    out.write_text(
        json.dumps(records, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def main(out, scout=False):
    job = next(j for j in json.loads(TRAIN.read_text())["jobs"] if j["id"] == 3142)
    p = np.array(job["parameters"])
    source = sample(np.array(job["initial"]), p, 200.0, 6000.0, 0.01)
    disp = source[:, 8:]
    row = int(np.argmin(np.r_[disp.mean(axis=0), -disp.sum(axis=1).mean()]))
    trees = {
        level: templates(section(source, p, level, 0.1)[1], row) for level in LEVELS
    }
    independent = np.array(json.loads(LARGE.read_text())["jobs"][0]["initial"])
    if scout:
        return screen(out, p, source, independent, trees[0.0])
    findings = []
    for b in (1.5075, 1.515, 1.52):
        p[2] = b
        checks = []
        for x0 in (source[-1], independent):
            states = sample(x0, p, 1000.0, 6000.0, 0.01)
            seq = labels(distances(trees[0.0], section(states, p, 0.0, 0.1)[1]), 2)
            assert len(seq) > 20
            m = next(k for k in range(1, 12) if np.all(seq[k:] == seq[:-k]))
            checks.append(
                {
                    "start": x0.tolist(),
                    "changes": len(seq) - 1,
                    "repeating_word": seq[:m].tolist(),
                }
            )
        sensitivities = []
        for level in LEVELS:
            d = distances(trees[level], section(states, p, level, 0.1)[1])
            for ratio in (1.5, 2, 2.5):
                seq = labels(d, ratio)
                word = None
                if len(seq) > 20:
                    for ell in range(1, 21):
                        if np.array_equal(seq[ell:], seq[:-ell]):
                            word = seq[:ell].tolist()
                            break
                sensitivities.append(
                    {
                        "level": level,
                        "ratio": ratio,
                        "labels": sorted(set(seq.tolist())),
                        "word": word,
                    }
                )
        z, period, diagnostics = get_cycle(states, p)
        grid = np.linspace(0, period, m + 1)
        verify = solve_ivp(
            lambda t, x: source_force_constraint_field(x, p, 4),
            (0, period),
            z,
            t_eval=grid,
            method="DOP853",
            rtol=2e-12,
            atol=2e-14,
            max_step=0.08,
        )
        assert verify.success
        errors = sorted(
            (
                max(
                    np.max(np.abs(act(x, perm) - y))
                    for x, y in zip(verify.y.T[:-1], verify.y.T[1:])
                ),
                perm,
            )
            for perm in permutations(range(4))
        )
        diagnostics.update(
            source_return_inf=float(np.max(np.abs(verify.y[:, -1] - z))),
            shift_permutation=list(errors[0][1]),
            shift_error_inf=float(errors[0][0]),
            second_permutation_error_inf=float(errors[1][0]),
        )
        findings.append(
            dict(
                b=b,
                period=period,
                state=z.tolist(),
                checks=checks,
                section_controls=sensitivities,
                **diagnostics,
            )
        )
        print(b, m, period, diagnostics, flush=True)
    out.write_text(
        json.dumps(findings, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--scout",
        action="store_true",
        help="Repeat the two parameter sweeps used to select the three regular itineraries.",
    )
    args = parser.parse_args()
    main(args.output, args.scout)

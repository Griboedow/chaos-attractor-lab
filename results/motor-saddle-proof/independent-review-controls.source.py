"""Independent bounded review controls; does not certify an orbit.

Execute with .venv-proof/Scripts/python.exe -B from the repository root.
The only output is a new adjacent JSON file. Existing evidence is preserved.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research"))
import flint
from flint import arb, arb_mat, ctx
import validate_motor_saddle_cycle as saddle
from validate_motor_cycle import ball_record, read_ball


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    output = Path(__file__).with_name("independent-review-controls.json")
    if output.exists():
        raise FileExistsError("Preserve prior review controls")
    names = ("validate_motor_saddle_cycle.py", "validate_motor_paired_cycle.py",
             "validate_motor_cycle.py", "check_motor_saddle_proof.py")
    sources = {name: (ROOT / "research" / name).read_bytes() for name in names}
    ctx.prec = 256
    checks = []
    for seed in (1, 7):
        state = [arb(((j * 7 + seed * 11) % 37) - 18) / 29
                 for j in range(saddle.DIM)]
        ode = saddle.SaddleODE(arb("817/250"), arb("237/50"),
                               with_transverse=True, order=8)
        initial = ode.initial(state)
        fast = ode.jets(initial, arb("1/7"), 8)
        ode.generic_jets = True
        generic = ode.jets(initial, arb("1/7"), 8)
        difference = max((a - b).abs_upper()
                         for u, v in zip(fast, generic) for a, b in zip(u, v))
        overlap = all(a.overlaps(b)
                      for u, v in zip(fast, generic) for a, b in zip(u, v))
        checks.append(dict(name="nonconstant matrix/scalar Taylor comparison",
                           seed=seed, components=sum(map(len, fast)),
                           interval_overlap=overlap,
                           difference=ball_record(difference),
                           passed=bool(overlap and difference < arb("1e-70"))))

    # Exercise production shooting assembly using an exactly soluble affine
    # field with A^2=0 and A*b=0. All numbers below are exact dyadic rationals;
    # centered differences are exact because the residual is affine in each
    # variable separately. This checks layout, phase removal and all 68
    # duration/state columns, independently of the motor ODE integrator.
    d = saddle.DIM
    a = arb_mat(d, d)
    b = [arb(0)] * d
    for i in range(8):
        a[i, 8 + i] = (-1)**i * (i + 1)
        b[i] = i - 3
    eye = saddle.identity(d)

    class AffineODE:
        def __init__(self, *args, with_transverse=False, **kwargs):
            self.offsets = {"core": (d, d)}
            if with_transverse:
                self.offsets.update({"2": (d + d*d, 3),
                                     "5": (d + d*d + 9, 3)})

        def value(self, state, duration):
            v = a * arb_mat([[x] for x in state])
            return [duration * (v[i, 0] + b[i]) for i in range(d)]

    def affine_flow(ode, state, duration, *unused):
        v = ode.value(state, arb(1))
        end = [x + duration * y for x, y in zip(state, v)]
        matrix = eye + duration * a
        return end + [matrix[i, j] for i in range(d) for j in range(d)]

    original_ode, original_flow = saddle.SaddleODE, saddle.flow
    try:
        saddle.SaddleODE, saddle.flow = AffineODE, affine_flow
        q = [arb(j - 31) / 8 for j in range(67)] + [arb("1/8")]
        residual, jac = saddle.shooting(q, arb(3), arb(4), 1, 2, segments=4)
        epsilon = arb("1/16")
        exact_columns = []
        for column in range(68):
            plus, minus = list(q), list(q)
            plus[column] += epsilon
            minus[column] -= epsilon
            fp = saddle.shooting(plus, arb(3), arb(4), 1, 2, segments=4)[0]
            fm = saddle.shooting(minus, arb(3), arb(4), 1, 2, segments=4)[0]
            exact_columns.extend(jac[row, column] == (fp[row, 0] - fm[row, 0]) / (2 * epsilon)
                                 for row in range(68))
        states = [q[:12] + [arb(0)] + q[12:16]]
        states += [q[16 + 17*k:33 + 17*k] for k in range(3)]
        p = saddle.permutation_matrix()
        expected = []
        for k, state in enumerate(states):
            endpoint = affine_flow(AffineODE(), state, q[-1])[:d]
            target = states[k + 1] if k < 3 else [v[0] for v in (p * arb_mat([[x] for x in states[0]])).tolist()]
            expected += [x - y for x, y in zip(endpoint, target)]
        exact_residual = [residual[i, 0] == expected[i] for i in range(68)]
        checks.append(dict(name="four-segment residual and 68-column Jacobian layout",
                           exact_derivative_entries=len(exact_columns),
                           exact_residual_entries=len(exact_residual),
                           passed=bool(all(exact_columns) and all(exact_residual))))

        # Noncommuting segment matrices distinguish forward composition from
        # the common reversed-product mistake. Spectrum is called only to
        # inspect assembly; its Floquet pass flag is intentionally irrelevant.
        core, pairs = [], {"2": [], "5": []}
        for k in range(4):
            m = saddle.identity(d)
            m[k % 3, (k + 1) % 3] = arb(k + 1) / 8
            core.append(m)
            for name in pairs:
                n = saddle.identity(3)
                n[k % 3, (k + 1) % 3] = arb(k + int(name)) / 16
                pairs[name].append(n)
        cursor = [0]

        def matrix_flow(ode, state, duration, *unused):
            k = cursor[0]
            cursor[0] += 1
            matrices = [core[k], pairs["2"][k], pairs["5"][k]]
            return [arb(0)] * d + [m[i, j] for m in matrices
                for i in range(m.nrows()) for j in range(m.ncols())]

        saddle.flow = matrix_flow
        result = saddle.spectrum(q, arb(3), arb(4), 1, 2, [], segments=4)
        half = core[3] * core[2] * core[1] * core[0]
        expected_core = (p * half)**2
        observed_core = arb_mat([[read_ball(x) for x in row]
                                for row in result["quotient_full_monodromy"]])
        pair_checks = []
        for name, matrices in pairs.items():
            expected_pair = (matrices[3] * matrices[2] * matrices[1] * matrices[0])**2
            observed_pair = arb_mat([[read_ball(x) for x in row]
                                    for row in result["transverse"][name]["monodromy"]])
            pair_checks.append(observed_pair == expected_pair)
        reverse = (p * core[0] * core[1] * core[2] * core[3])**2
        checks.append(dict(name="noncommuting core and pair segment composition",
                           reverse_order_distinguished=bool(reverse != expected_core),
                           passed=bool(observed_core == expected_core and all(pair_checks)
                                       and reverse != expected_core)))
    finally:
        saddle.SaddleODE, saddle.flow = original_ode, original_flow

    attempt_path = ROOT / "results/motor-saddle-proof/attempt-01.json"
    attempt = json.loads(attempt_path.read_text())
    snapshot_matches = {name: digest(attempt_path.with_name(attempt_path.stem + "." + name).read_bytes()) == sha
                        for name, sha in attempt["source_sha256"].items()}
    seed = json.loads((ROOT / "results/motor-saddle-proof/multiple-shooting-seed.json").read_text())
    source_matches = digest((ROOT / seed["source"]).read_bytes()) == seed["source_sha256"]
    checks.append(dict(name="failed-attempt source snapshots and multiple-shooting seed provenance",
                       snapshots=snapshot_matches, seed_npz_matches=source_matches,
                       passed=bool(all(snapshot_matches.values()) and source_matches)))
    if any((ROOT / "research" / name).read_bytes() != data for name, data in sources.items()):
        raise RuntimeError("A reviewed source changed while controls were running")
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  status="passed" if all(x["passed"] for x in checks) else "failed",
                  precision_bits=ctx.prec, python=sys.version, python_flint=flint.__version__,
                  source_sha256={name: digest(data) for name, data in sources.items()},
                  helper_sha256=digest(Path(__file__).read_bytes()), checks=checks,
                  limitation="Bounded implementation controls and source review; not an orbit existence or Floquet certificate.")
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(record["status"], output.name)
    if record["status"] != "passed":
        raise SystemExit(2)


if __name__ == "__main__":
    main()

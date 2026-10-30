"""Independent short step control from the exact inherited CPU QR endpoint."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.integrate import solve_ivp


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, default=Path(__file__).parent / 'run-01')
    parser.add_argument('--output', type=Path, default=Path(__file__).parent / 'strict-control.json')
    args = parser.parse_args()
    run, out = args.run.resolve(), args.output.resolve()
    if out.exists():
        raise FileExistsError('Preserve existing short-control evidence')
    spec = importlib.util.spec_from_file_location('frozen_motor_rows', run / 'motor_rows.py')
    field = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = field
    spec.loader.exec_module(field)
    with np.load(run / 'input.npz') as data:
        p = data['parameters'].copy()
        initial = data['final_state'].copy()
        initial_basis = data['final_basis'].copy()
    with np.load(run / 'block-01.npz') as data:
        recorded = data['states'][:50].copy()
    D, N, dt = 23, 8, .2
    theta, constant = 7 / 8, -16 + 7 * (p[1] - p[0])

    def geometry(x):
        c, s, positions = field.unpack(x, N)
        v = -p[0] * positions + p[1] * (s - s.mean())
        return theta * (p[1] * c.sum() + p[0] * np.dot(positions, positions) / 2), np.dot(v, v)

    def augmented(t, state):
        f, j = field.field_jac(state[:D], p, N)
        return np.r_[f, (j @ state[D:D + D * D].reshape(D, D)).ravel(), np.trace(j), geometry(state[:D])[1]]

    def integrate(step, rtol, atol):
        x, basis = initial.copy(), initial_basis.copy()
        logs, trace, speed = np.zeros(D), 0., 0.
        states = []
        started = time.perf_counter()
        for _ in range(50):
            solution = solve_ivp(augmented, (0, dt), np.r_[x, basis.ravel(), 0., 0.],
                                 method='DOP853', rtol=rtol, atol=atol, max_step=step)
            assert solution.success
            x = solution.y[:D, -1]
            basis, triangular = np.linalg.qr(solution.y[D:D + D * D, -1].reshape(D, D))
            logs += np.log(abs(np.diag(triangular)))
            trace += solution.y[-2, -1]
            speed += solution.y[-1, -1]
            states.append(x.copy())
        residual = trace + geometry(x)[0] - geometry(initial)[0] - constant * 10 + theta * speed
        return dict(states=np.asarray(states), basis=basis, logs=logs, trace=trace, speed=speed,
                    identity_residual=float(residual), wall_seconds=time.perf_counter() - started,
                    settings=dict(duration_tau=10., max_step=step, rtol=rtol, atol=atol, qr_tau=dt))

    baseline = integrate(.03, 1e-11, 1e-13)
    strict = integrate(.015, 2e-13, 2e-15)
    errors = dict(baseline_recorded_states=float(np.max(abs(baseline['states'] - recorded))),
                  strict_recorded_states=float(np.max(abs(strict['states'] - recorded))),
                  final_basis=float(np.max(abs(strict['basis'] - baseline['basis']))),
                  log_sums=float(np.max(abs(strict['logs'] - baseline['logs']))),
                  trace=abs(float(strict['trace'] - baseline['trace'])),
                  speed_integral=abs(float(strict['speed'] - baseline['speed'])))
    controls = dict(exact_endpoint_basis_restored=True,
                    saved_trajectory_reproduced=errors['baseline_recorded_states'] < 1e-10,
                    step_tolerance_agreement=max(errors.values()) < 1e-9,
                    baseline_weighted_identity=abs(baseline['identity_residual']) < 1e-9,
                    strict_weighted_identity=abs(strict['identity_residual']) < 1e-9)
    assert all(controls.values()), (controls, errors)
    for result in (baseline, strict):
        for name in ('states', 'basis', 'logs'):
            result[name] = result[name].tolist()
        result['trace'], result['speed'] = float(result['trace']), float(result['speed'])
    record = dict(status='passed', created_utc=datetime.now(timezone.utc).isoformat(), controls=controls,
                  errors=errors, baseline=baseline, strict=strict,
                  source_sha256=sha(Path(__file__)),
                  input_hashes={name: sha(run / name) for name in ('input.npz', 'block-01.npz', 'motor_rows.py')},
                  scope='Ten original-tau units from the frozen 4000-observed endpoint, with all 23 QR directions. This local step/tolerance check does not establish a global trajectory error bound over 20000 tau.')
    out.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=record['status'], errors=errors)))


if __name__ == '__main__':
    main()

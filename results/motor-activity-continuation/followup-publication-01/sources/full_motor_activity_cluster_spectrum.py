"""Complete conditional tangent spectrum of the archived lane-46 cluster state.

All 14 quotient directions and three independent three-dimensional contrast
blocks are retained. The declared restriction is not a full-space orbit probe.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag
import motor_rows
import motor_postcritical_cycles as cycles

ROOT = Path(__file__).resolve().parents[1]
GROUPS = [[0], [1, 2], [3, 6], [4, 7], [5]]
WEIGHTS = np.array([1, 2, 2, 2, 1]) / 8
PAIRED = [1, 2, 3]
D = 14


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')


def contrast(state, parameters, j):
    rho, r, _ = parameters
    c, s, x = cycles.unpack(state, WEIGHTS)
    v = -rho * x + r * (s - WEIGHTS @ s)
    return np.array([[-1, -v[j] - r * s[j], rho * s[j]],
                     [v[j], -1 + r * c[j], -rho * c[j]], [0, r, -rho]])


def embeddings():
    quotient = np.column_stack([cycles.lift(v, GROUPS) for v in np.eye(D)])
    blocks = [quotient]
    for j in PAIRED:
        a, b = GROUPS[j]
        matrix = np.zeros((23, 3))
        matrix[a, 0], matrix[b, 0] = 1, -1
        matrix[8 + a, 1], matrix[8 + b, 1] = 1, -1
        if a < 7:
            matrix[16 + a, 2] = 1
        if b < 7:
            matrix[16 + b, 2] = -1
        blocks.append(matrix)
    return np.column_stack(blocks)


def geometry(state, parameters):
    c, s, x = cycles.unpack(state, WEIGHTS)
    v = -parameters[0] * x + parameters[1] * (s - WEIGHTS @ s)
    phi = 7 / 8 * (parameters[1] * (8 * WEIGHTS @ c) + parameters[0] * (8 * WEIGHTS @ (x * x)) / 2)
    return float(phi), float(8 * WEIGHTS @ (v * v))


def decomposition_control(state, parameters, transform):
    full_state = cycles.lift(state, GROUPS)
    f, jq = cycles.field_jac(state, parameters, WEIGHTS)
    ff, jf = motor_rows.field_jac(full_state, parameters, 8)
    blocks = [jq] + [contrast(state, parameters, j) for j in PAIRED]
    errors = dict(field=float(np.max(abs(transform[:, :D] @ f - ff))),
                  all_23_directions=float(np.max(abs(jf @ transform - transform @ block_diag(*blocks)))),
                  trace=float(abs(np.trace(jf) - sum(np.trace(m) for m in blocks))),
                  independent_force=float(np.max(abs(ff - motor_rows.source_force_constraint_field(full_state, parameters, 8)))))
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--observed', type=float, default=4000.)
    parser.add_argument('--wall', type=float, default=1200.)
    parser.add_argument('--max-step', type=float, default=.03)
    parser.add_argument('--rtol', type=float, default=1e-11)
    parser.add_argument('--atol', type=float, default=1e-13)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve previous numerical evidence')
    if args.observed <= 0 or args.observed % 500 or min(args.wall, args.max_step, args.rtol, args.atol) <= 0:
        raise ValueError('Positive settings and observation a multiple of 500 required')
    source = ROOT / 'results/motor-rows/activity-path/N8-lifetime-history.npz'
    parent_analysis = ROOT / 'results/motor-activity-continuation/cluster-controls/lane-46/analysis.json'
    expected = json.loads(parent_analysis.read_text())
    with np.load(source) as data:
        lane = int(np.flatnonzero(data['source_lane_indices'] == 46)[0])
        original = data['final_states'][lane].copy()
        parameters = data['parameters'][lane].copy()
    state = cycles.restrict(original, GROUPS)
    projection_error = float(np.max(abs(cycles.lift(state, GROUPS) - original)))
    assert GROUPS == expected['groups'] and np.array_equal(WEIGHTS, expected['weights'])
    assert np.array_equal(parameters, expected['parameters']) and sha(source) == expected['source_sha256']
    assert projection_error == expected['initial_projection_error']
    transform = embeddings()
    assert np.linalg.matrix_rank(transform) == 23
    controls = [decomposition_control(state, parameters, transform)]
    rng = np.random.default_rng(2026103046)
    controls += [decomposition_control(rng.normal(size=D), parameters, transform) for _ in range(12)]
    assert max(value for item in controls for value in item.values()) < 1e-10
    restriction = cycles.audit(GROUPS, parameters)
    out = args.output.resolve()
    out.mkdir(parents=True)
    shutil.copyfile(parent_analysis, out / 'parent-2000-analysis.json')
    source_hashes = {}
    for path in (Path(__file__), Path(cycles.__file__), Path(motor_rows.__file__)):
        shutil.copyfile(path, out / path.name)
        source_hashes[path.name] = sha(out / path.name)
    np.savez_compressed(out / 'seed.npz', original_full_state=original, restricted_state=state,
                        parameters=parameters, weights=WEIGHTS, full_tangent_transform=transform)
    record = dict(status='running', started_utc=datetime.now(timezone.utc).isoformat(), source_lane=46,
                  source_path=source.relative_to(ROOT).as_posix(), source_sha256=sha(source),
                  parent_analysis_sha256=sha(parent_analysis), seed_sha256=sha(out / 'seed.npz'),
                  source_hashes=source_hashes, groups=GROUPS, weights=WEIGHTS.tolist(), parameters=parameters.tolist(),
                  initial_projection_error=projection_error, restriction_controls=restriction,
                  complete_decomposition_controls=controls,
                  transform_determinant=float(np.linalg.det(transform)), transform_condition=float(np.linalg.cond(transform)),
                  settings=dict(warm_tau=200., observed_tau_requested=args.observed, block_tau=500., qr_tau=.5,
                                max_step=args.max_step, rtol=args.rtol, atol=args.atol, wall_seconds=args.wall),
                  runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, blas_threads=1),
                  tangent_blocks=[dict(kind='quotient', dimension=14)] + [dict(kind='pair_contrast', rows=GROUPS[j], dimension=3) for j in PAIRED],
                  max_local_trace_errors=[0., 0., 0., 0.], max_local_flow_error=0., max_condition=0.,
                  windows=[], interpretation='A declared weighted invariant restriction with a measured tiny initial projection. All 23 tangent directions are represented by a constant invertible decomposition. Finite-time QR values use separate block Euclidean norms; their asymptotic union is the full tangent spectrum if limits exist.',
                  limitations=['No full-space physical perturbation in this run.', 'Finite-time estimates do not prove asymptotic chaos, uniform transverse stability, or attraction.'])
    sizes = [D, 3, 3, 3]
    offsets = np.cumsum([D] + [size * size for size in sizes])
    matrices = [np.eye(size) for size in sizes]
    logs = [np.zeros(size) for size in sizes]
    trace = np.zeros(4)
    times, states, checkpoints = [], [], []
    observed = elapsed = flow_log = speed_integral = 0.
    phi_start = flow_initial = None
    constant = -16 + 7 * (parameters[1] - parameters[0])
    began = time.perf_counter()

    def augmented(t, z):
        x = z[:D]
        f, jq = cycles.field_jac(x, parameters, WEIGHTS)
        jacobians = [jq] + [contrast(x, parameters, j) for j in PAIRED]
        parts = [f] + [(j @ z[a:b].reshape(size, size)).ravel()
                       for j, a, b, size in zip(jacobians, offsets[:-1], offsets[1:], sizes)]
        return np.r_[np.concatenate(parts), [np.trace(j) for j in jacobians], geometry(x, parameters)[1]]

    # Check the complete coordinate decomposition at the level of a finite-time
    # tangent map, independently integrating the original full 23D equation.
    initial_augmented = np.r_[state, np.eye(D).ravel(), np.tile(np.eye(3).ravel(), 3), np.zeros(5)]
    short = solve_ivp(augmented, (0, 1.), initial_augmented, method='DOP853', rtol=2e-13, atol=2e-15, max_step=.015)
    def full_augmented(t, z):
        f, j = motor_rows.field_jac(z[:23], parameters, 8)
        return np.r_[f, (j @ z[23:].reshape(23, 23)).ravel()]
    full = solve_ivp(full_augmented, (0, 1.), np.r_[cycles.lift(state, GROUPS), transform.ravel()],
                     method='DOP853', rtol=2e-13, atol=2e-15, max_step=.015)
    assert short.success and full.success
    end = short.y[:, -1]
    block_map = block_diag(*[end[a:b].reshape(size, size) for a, b, size in zip(offsets[:-1], offsets[1:], sizes)])
    short_errors = dict(state=float(np.max(abs(cycles.lift(end[:D], GROUPS) - full.y[:23, -1]))),
                        all_tangent_columns=float(np.max(abs(transform @ block_map - full.y[23:, -1].reshape(23, 23)))))
    assert max(short_errors.values()) < 1e-10
    record['independent_full_23D_short_map'] = dict(duration_tau=1., errors=short_errors, passed=True)

    def save():
        if observed:
            rates = [part / observed for part in logs]
            record.update(quotient_spectrum=rates[0].tolist(), contrast_spectra=[r.tolist() for r in rates[1:]],
                          sorted_union_spectrum=np.sort(np.concatenate(rates))[::-1].tolist(),
                          integrated_block_traces=trace.tolist(),
                          block_trace_errors_per_tau=[float(abs(part.sum() - t) / observed) for part, t in zip(logs, trace)],
                          weighted_volume_residual=float(trace.sum() + geometry(state, parameters)[0] - phi_start - constant * observed + 7 / 8 * speed_integral),
                          flow_exponent=flow_log / observed,
                          telescopic_flow_exponent=float(np.log(np.linalg.norm(cycles.field_jac(state, parameters, WEIGHTS)[0]) / flow_initial) / observed))
        record.update(elapsed_tau=elapsed, observed_tau=observed, wall_seconds=time.perf_counter() - began,
                      final_decomposition_control=decomposition_control(state, parameters, transform))
        np.savez_compressed(out / 'history.npz', time=np.asarray(times), states=np.asarray(states),
                            parameters=parameters, weights=WEIGHTS, final_state=state,
                            final_quotient_basis=matrices[0], final_contrast_bases=np.asarray(matrices[1:]),
                            log_sums=np.concatenate(logs), integrated_block_traces=trace,
                            integrated_speed_squared=speed_integral, normalized_flow_log=flow_log,
                            checkpoints=np.asarray(checkpoints))
        record['history_sha256'] = sha(out / 'history.npz')
        write(out / 'analysis.json', record)

    try:
        for block, length in enumerate([200.] + [500.] * int(args.observed / 500)):
            block_logs = [np.zeros(size) for size in sizes]
            block_trace = np.zeros(4)
            spent = block_flow = block_speed = 0.
            for _ in range(int(length * 2)):
                if time.perf_counter() - began >= args.wall:
                    record['status'] = 'wall budget incomplete'
                    break
                f0 = cycles.field_jac(state, parameters, WEIGHTS)[0]
                norm0 = np.linalg.norm(f0)
                z0 = np.r_[state, np.concatenate([m.ravel() for m in matrices]), np.zeros(5)]
                solution = solve_ivp(augmented, (0, .5), z0, method='DOP853', rtol=args.rtol,
                                     atol=args.atol, max_step=args.max_step)
                if not solution.success:
                    raise ArithmeticError(solution.message)
                end = solution.y[:, -1]
                raw = [end[a:b].reshape(size, size) for a, b, size in zip(offsets[:-1], offsets[1:], sizes)]
                f1 = cycles.field_jac(end[:D], parameters, WEIGHTS)[0]
                transported = raw[0] @ (matrices[0].T @ (f0 / norm0))
                flow_error = float(np.linalg.norm(transported * norm0 - f1) / (1 + np.linalg.norm(f1)))
                record['max_local_flow_error'] = max(record['max_local_flow_error'], flow_error)
                assert flow_error < 1e-8
                for k, matrix in enumerate(raw):
                    q, r = np.linalg.qr(matrix)
                    growth = np.log(abs(np.diag(r)))
                    error = float(abs(growth.sum() - end[offsets[-1] + k]))
                    condition = float(np.linalg.cond(matrix))
                    assert error < 1e-8 and condition < 1e10
                    record['max_local_trace_errors'][k] = max(record['max_local_trace_errors'][k], error)
                    record['max_condition'] = max(record['max_condition'], condition)
                    block_logs[k] += growth
                    matrices[k] = q
                block_trace += end[offsets[-1]:offsets[-1] + 4]
                block_flow += float(np.log(np.linalg.norm(transported)))
                block_speed += float(end[-1])
                state = end[:D]
                elapsed += .5
                spent += .5
                times.append(elapsed)
                states.append(state.copy())
            if block == 0:
                phi_start = geometry(state, parameters)[0]
                flow_initial = float(np.linalg.norm(cycles.field_jac(state, parameters, WEIGHTS)[0]))
            elif spent:
                observed += spent
                for total, part in zip(logs, block_logs):
                    total += part
                trace += block_trace
                flow_log += block_flow
                speed_integral += block_speed
                record['windows'].append(dict(end_tau=elapsed, duration=spent, complete=spent == length,
                                               quotient_spectrum=(block_logs[0] / spent).tolist(),
                                               contrast_spectra=[(part / spent).tolist() for part in block_logs[1:]],
                                               integrated_traces=block_trace.tolist()))
            checkpoints.append(np.r_[elapsed, state, np.concatenate([m.ravel() for m in matrices]),
                                      np.concatenate(logs), trace, speed_integral, flow_log])
            save()
            print(json.dumps(dict(elapsed_tau=elapsed, observed_tau=observed,
                                  quotient=record.get('quotient_spectrum'), contrasts=record.get('contrast_spectra'),
                                  wall_seconds=round(record['wall_seconds'], 2))), flush=True)
            if record['status'] != 'running':
                break
        if record['status'] == 'running':
            record['status'] = 'complete'
    except Exception as error:
        record.update(status='error', error=repr(error))
        raise
    finally:
        save()
    if record['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    main()

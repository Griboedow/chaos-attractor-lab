"""Classify both excluded tuple-15 endpoints without discarding a transient."""
import argparse
from contextlib import redirect_stdout
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.integrate import solve_ivp
import motor_rows
import motor_postcritical_cycles as cycles

ROOT = Path(__file__).resolve().parents[1]
D = 23


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')


def rows(state):
    return np.column_stack((state[:8], state[8:16], np.r_[state[16:], -state[16:].sum()]))


def groups(state):
    points = rows(state)
    todo, result = set(range(8)), []
    while todo:
        first = min(todo)
        group = [j for j in sorted(todo) if np.linalg.norm(points[j] - points[first]) < 1e-10]
        result.append(group)
        todo -= set(group)
    return result


def lags(times, states):
    result = []
    points = np.asarray(states)
    for lag in range(1, min(16, len(points) - 1) + 1):
        errors = np.max(abs(points[lag:] - points[:-lag]), axis=1)
        tail = errors[-min(12, len(errors)):]
        result.append(dict(lag=lag, period=float(times[-1] - times[-1 - lag]),
                           final_error=float(errors[-1]), tail_max_error=float(tail.max()),
                           tail_rms_error=float(np.sqrt(np.mean(tail * tail))), tail_count=len(tail)))
    return result


def endpoint_control(state, p):
    f, j = motor_rows.field_jac(state, p, 8)
    complex_j = np.column_stack([np.imag(motor_rows.rhs(0, state.astype(complex) + 1e-27j * v, p, 8)) / 1e-27 for v in np.eye(D)])
    return dict(jacobian=float(np.max(abs(j - complex_j))),
                independent_force=float(np.max(abs(f - motor_rows.source_force_constraint_field(state, p, 8)))),
                trace=float(abs(np.trace(j) - (-16 - 7 * p[0] + 7 / 8 * p[1] * state[:8].sum()))))


def run(lane, out, duration, wall):
    source = ROOT / 'results/motor-rows/activity-path/N8-refined.npz'
    scout_path = source.with_name('N8.npz')
    with np.load(source) as data, np.load(scout_path) as scout:
        index = int(np.flatnonzero(data['source_lane_indices'] == lane)[0])
        initial, tangent, p = [data[key][index].copy() for key in ('final_states', 'final_tangents', 'parameters')]
        inherited_segment = float(data['elapsed_original'][index])
        scout_elapsed = float(scout['elapsed_original'][lane])
        assert np.array_equal(data['initial_states'][index], scout['final_states'][lane])
        assert int(data['original_tuple_id'][index]) == 15
        original_windows = data['lle_windows'][index].tolist()
        start = int(data['original_start_id'][index])
    out.mkdir(parents=True)
    cumulative_start = scout_elapsed + inherited_segment
    normalized = tangent / np.linalg.norm(tangent)
    u = np.eye(D)[0] - normalized
    basis = np.eye(D) if np.linalg.norm(u) < 1e-14 else np.eye(D) - 2 * np.outer(u, u) / np.dot(u, u)
    assert np.max(abs(basis[:, 0] - normalized)) < 1e-12
    assert np.max(abs(basis.T @ basis - np.eye(D))) < 1e-12
    state, logs, trace, flow_log = initial.copy(), np.zeros(D), 0., 0.
    initial_flow = np.linalg.norm(motor_rows.rhs(0, state, p, 8))
    record = dict(status='running', started_utc=datetime.now(timezone.utc).isoformat(), source_lane=lane,
                  original_tuple_id=15, original_start_id=start, parameters=p.tolist(), initial_state=initial.tolist(),
                  input_source_sha256=sha(source), scout_source_sha256=sha(scout_path),
                  source_path=source.relative_to(ROOT).as_posix(), scout_path=scout_path.relative_to(ROOT).as_posix(),
                  original_refined_windows=original_windows, refined_segment_elapsed=inherited_segment,
                  scout_elapsed=scout_elapsed, cumulative_source_elapsed=cumulative_start,
                  refined_initial_state_equals_scout_endpoint=True, initial_groups=groups(initial),
                  initial_endpoint_control=endpoint_control(initial, p), source_hashes={}, windows=[],
                  settings=dict(duration_tau=duration, warm_tau=0., qr_tau=.5, max_step=.025, rtol=2e-12, atol=2e-14,
                                wall_seconds=wall, state_projection=False),
                  tangent_initialization='Carry the normalized archived leading tangent as the first vector; complete it by a deterministic Householder orthonormal basis. Other directions and all reported sums start at this diagnostic endpoint; no inherited full basis exists.',
                  initial_carried_tangent_error=float(np.max(abs(basis[:, 0] - normalized))),
                  max_local_trace_error=0., max_local_flow_error=0., max_condition=0.)
    assert max(record['initial_endpoint_control'].values()) < 1e-10
    for path in (Path(__file__), Path(motor_rows.__file__), Path(cycles.__file__)):
        shutil.copyfile(path, out / path.name)
        record['source_hashes'][path.name] = sha(out / path.name)
    np.savez_compressed(out / 'seed.npz', initial_state=initial, archived_tangent=tangent,
                        initial_basis=basis, parameters=p)
    record['seed_sha256'] = sha(out / 'seed.npz')
    times, states, section_times, section_states = [0.], [initial.copy()], [], []
    elapsed = 0.
    began = time.perf_counter()

    def augmented(t, z):
        f, j = motor_rows.field_jac(z[:D], p, 8)
        return np.r_[f, (j @ z[D:D + D * D].reshape(D, D)).ravel(), np.trace(j)]

    def section(t, z):
        return z[16]
    section.direction = 1

    def save():
        record.update(elapsed_added_tau=elapsed, cumulative_elapsed_tau=cumulative_start + elapsed,
                      wall_seconds=time.perf_counter() - began, final_groups=groups(state),
                      final_endpoint_control=endpoint_control(state, p), section_count=len(section_times),
                      return_lags=lags(section_times, section_states))
        if elapsed:
            record.update(spectrum=(logs / elapsed).tolist(), trace_error_per_tau=float(abs(logs.sum() - trace) / elapsed),
                          flow_exponent=flow_log / elapsed,
                          telescopic_flow_exponent=float(np.log(np.linalg.norm(motor_rows.rhs(0, state, p, 8)) / initial_flow) / elapsed))
        np.savez_compressed(out / 'sample.npz', parameters=p, time=np.array(times), states=np.array(states),
                            section_time=np.array(section_times), section_states=np.array(section_states).reshape(-1, D),
                            initial_state=initial, final_state=state, final_basis=basis,
                            log_sums=logs, integrated_trace=trace, normalized_flow_log=flow_log)
        record['sample_sha256'] = sha(out / 'sample.npz')
        write(out / 'analysis.json', record)

    try:
        for block in range(int(duration / 100)):
            window_logs, spent = np.zeros(D), 0.
            for _ in range(200):
                if time.perf_counter() - began >= wall:
                    record['status'] = 'wall budget incomplete'
                    break
                f0 = motor_rows.rhs(0, state, p, 8)
                norm = np.linalg.norm(f0)
                sol = solve_ivp(augmented, (0, .5), np.r_[state, basis.ravel(), 0.], method='DOP853',
                                 rtol=2e-12, atol=2e-14, max_step=.025, events=section)
                if not sol.success:
                    raise ArithmeticError(sol.message)
                end = sol.y[:D, -1]
                matrix = sol.y[D:D + D * D, -1].reshape(D, D)
                new_basis, triangular = np.linalg.qr(matrix)
                growth = np.log(abs(np.diag(triangular)))
                delta_trace = float(sol.y[-1, -1])
                error = float(abs(growth.sum() - delta_trace))
                propagated = matrix @ (basis.T @ (f0 / norm))
                f1 = motor_rows.rhs(0, end, p, 8)
                flow_error = float(np.linalg.norm(propagated * norm - f1) / (1 + np.linalg.norm(f1)))
                condition = float(np.linalg.cond(matrix))
                assert error < 1e-8 and flow_error < 1e-8 and condition < 1e10
                record['max_local_trace_error'] = max(record['max_local_trace_error'], error)
                record['max_local_flow_error'] = max(record['max_local_flow_error'], flow_error)
                record['max_condition'] = max(record['max_condition'], condition)
                for hit, z in zip(sol.t_events[0], sol.y_events[0]):
                    if hit > 1e-10:
                        section_times.append(elapsed + float(hit))
                        section_states.append(z[:D])
                state, basis = end, new_basis
                elapsed += .5
                spent += .5
                logs += growth
                window_logs += growth
                trace += delta_trace
                flow_log += float(np.log(np.linalg.norm(propagated)))
                times.append(elapsed)
                states.append(state.copy())
            if spent:
                record['windows'].append(dict(end_tau=elapsed, duration=spent, spectrum=(window_logs / spent).tolist(),
                                               return_lags=lags(section_times, section_states)))
            save()
            print(json.dumps(dict(lane=lane, elapsed=elapsed, spectrum=record.get('spectrum', [])[:6],
                                  best_return=min(record['return_lags'], key=lambda q:q['final_error']) if record['return_lags'] else None)), flush=True)
            if record['status'] != 'running':
                break
        if record['status'] == 'running':
            record['status'] = 'complete'
    except Exception as error:
        record.update(status='error', error=repr(error))
        raise
    finally:
        save()
    candidates = [v for v in record['return_lags'] if v['final_error'] < 1e-7 and v['tail_max_error'] < 1e-6]
    record['periodic_attempts'] = []
    if record['status'] == 'complete' and candidates:
        candidate = candidates[0]
        with (out / 'periodic-attempt.log').open('w', encoding='utf-8', newline='\n') as log_file:
            try:
                with redirect_stdout(log_file):
                    cycles.run(out / 'sample.npz', groups(section_states[-1]), candidate['lag'], out / 'periodic')
                result = json.loads((out / 'periodic/analysis.json').read_text())
                record['periodic_attempts'].append(dict(lag=candidate['lag'], status=result['status'],
                    period=result['period'], unstable_count=result['full_unstable_count_at_1e_minus_6'],
                    neutral_count=result['full_neutral_count_at_1e_minus_6']))
            except Exception as error:
                record['periodic_attempts'].append(dict(lag=candidate['lag'], status='failed', error=repr(error)))
    record['periodic_search_scope'] = 'All return lags1 through16 archived. Attempt the shortest close-return candidate, requiring final error<1e-7 and maximum over last up-to12 comparisons<1e-6. No failed lane discarded.'
    write(out / 'analysis.json', record)
    print(json.dumps(dict(lane=lane, status=record['status'], periodic_attempts=record['periodic_attempts'])), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--duration', type=float, default=500.)
    parser.add_argument('--wall-per-lane', type=float, default=600.)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve existing classification')
    if args.duration <= 0 or args.duration % 100:
        raise ValueError('Duration must be a positive multiple of100')
    args.output.mkdir(parents=True)
    for lane in (45, 47):
        run(lane, args.output / f'lane-{lane}', args.duration, args.wall_per_lane)


if __name__ == '__main__':
    main()

"""Bounded numerical continuation of the case-16 five-cluster periodic branch."""
import argparse
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
from scipy.linalg import block_diag
import motor_postcritical_cycles as cycles
import motor_rows

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')


def eig_record(values):
    return [dict(real=float(v.real), imag=float(v.imag), modulus=float(abs(v))) for v in values]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--wall', type=float, default=900.)
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError('Preserve previous continuation evidence')
    source = ROOT / 'results/motor-activity-continuation/probe-destinations/case-16/periodic'
    metadata = json.loads((source / 'analysis.json').read_text())
    with np.load(source / 'cycle.npz') as data:
        state = data['quotient_initial_state'].copy()
        p = data['parameters'].copy()
        period = float(data['time'][-1])
    groups = metadata['groups']
    weights = np.array(metadata['weights'])
    assert groups == [[0, 3], [1, 2], [4], [5, 6], [7]]
    assert np.array_equal(weights, np.array([2, 2, 1, 2, 1]) / 8)
    paired = [0, 1, 3]
    d = len(state)
    target = float(np.load(ROOT / 'results/motor-rows/activity-path/N8-refined.npz')['parameters'][
        int(np.flatnonzero(np.load(ROOT / 'results/motor-rows/activity-path/N8-refined.npz')['source_lane_indices'] == 46)[0]), 0])
    out.mkdir(parents=True)
    shutil.copyfile(source / 'analysis.json', out / 'input-analysis.json')
    shutil.copyfile(source / 'cycle.npz', out / 'input-cycle.npz')
    sources = {}
    for path in (Path(__file__), Path(cycles.__file__), Path(motor_rows.__file__)):
        shutil.copyfile(path, out / path.name)
        sources[path.name] = sha(out / path.name)
    record = dict(status='running', started_utc=datetime.now(timezone.utc).isoformat(), groups=groups,
                  weights=weights.tolist(), fixed_r_b=p[1:].tolist(), start_rho=float(p[0]), target_rho=target,
                  input_analysis_sha256=sha(out / 'input-analysis.json'), input_cycle_sha256=sha(out / 'input-cycle.npz'),
                  source_hashes=sources, settings=dict(maximum_step=.005, minimum_step=.0003125,
                  maximum_attempts=60, wall_seconds=args.wall, shooting_rtol=2e-11, shooting_max_step=.025,
                  spectral_rtol=2e-12, spectral_atol=2e-14, spectral_max_step=.02), attempts=[],
                  limitation='Floating-point continuation of one periodic branch, including unstable portions. No interval theorem, daughter-branch construction or connection to the lane46 chaotic trajectory is established.')
    transform_parts = [np.column_stack([cycles.lift(v, groups) for v in np.eye(d)])]
    for j in paired:
        a, b = groups[j]
        e = np.zeros((23, 3))
        e[a, 0], e[b, 0] = 1, -1
        e[8 + a, 1], e[8 + b, 1] = 1, -1
        if a < 7:
            e[16 + a, 2] = 1
        if b < 7:
            e[16 + b, 2] = -1
        transform_parts.append(e)
    transform = np.column_stack(transform_parts)
    assert np.linalg.matrix_rank(transform) == 23

    def contrast(x, parameters, j):
        c, s, positions = cycles.unpack(x, weights)
        rho, r, _ = parameters
        velocity = -rho * positions + r * (s - weights @ s)
        return np.array([[-1, -velocity[j] - r * s[j], rho * s[j]],
                         [velocity[j], -1 + r * c[j], -rho * c[j]], [0, r, -rho]])

    def spectral(x, t, parameters):
        def fun(tau, z):
            f, jq = cycles.field_jac(z[:d], parameters, weights)
            output = [f, (jq @ z[d:d + d*d].reshape(d, d)).ravel()]
            for k, j in enumerate(paired):
                start = d + d*d + 9*k
                output.append((contrast(z[:d], parameters, j) @ z[start:start+9].reshape(3, 3)).ravel())
            return np.concatenate(output)
        sol = solve_ivp(fun, (0, t), np.r_[x, np.eye(d).ravel(), np.tile(np.eye(3).ravel(), 3)],
                        method='DOP853', rtol=2e-12, atol=2e-14, max_step=.02)
        if not sol.success:
            raise ArithmeticError(sol.message)
        end = sol.y[:, -1]
        mq = end[d:d+d*d].reshape(d, d)
        blocks = end[d+d*d:].reshape(3, 3, 3)
        f, jq = cycles.field_jac(x, parameters, weights)
        _, jf = motor_rows.field_jac(cycles.lift(x, groups), parameters, 8)
        errors = dict(closure=float(np.max(abs(end[:d] - x))),
                      flow=float(np.max(abs(mq @ f - f))),
                      full_decomposition=float(np.max(abs(jf @ transform - transform @ block_diag(jq, *[contrast(x, parameters, j) for j in paired])))))
        assert max(errors.values()) < 1e-7
        values = np.linalg.eigvals(mq)
        flow_index = int(np.argmin(abs(values - 1)))
        nonflow = np.delete(values, flow_index)
        ordered = nonflow[np.argsort(-abs(nonflow))]
        transverse = [np.linalg.eigvals(m) for m in blocks]
        return mq, blocks, dict(errors=errors, quotient_flow_multiplier=eig_record([values[flow_index]])[0],
            quotient_nonflow_multipliers=eig_record(ordered), contrast_multipliers=[eig_record(v) for v in transverse],
            quotient_max_nonflow=float(abs(ordered[0])), contrast_maxima=[float(np.max(abs(v))) for v in transverse],
            quotient_unstable_count=int(np.sum(abs(nonflow) > 1 + 1e-6)),
            contrast_unstable_count=int(sum(np.sum(abs(v) > 1 + 1e-6) for v in transverse)))

    began = time.perf_counter()
    step, attempted_seed = .005, False
    try:
        for attempt in range(1, 61):
            if time.perf_counter() - began >= args.wall:
                record['status'] = 'wall budget incomplete'
                break
            rho = float(p[0]) if not attempted_seed else max(target, float(p[0]) - step)
            candidate_p = p.copy()
            candidate_p[0] = rho
            directory = out / f'attempt-{attempt:02d}'
            directory.mkdir()
            item = dict(index=attempt, rho=rho, step_from_previous=float(p[0] - rho),
                        input_period=period, input_state=state.tolist(), status='running')
            started = time.perf_counter()
            try:
                x, t, shot, history = cycles.shoot(state, period, candidate_p, weights)
                mq, blocks, spectrum = spectral(x, t, candidate_p)
                comparison = float(np.max(abs(mq - shot.y[d:, -1].reshape(d, d))))
                assert comparison < 1e-7
                item.update(status='complete', period=t, initial_state=x.tolist(), shooting_history=history,
                            spectrum=spectrum, strict_shooting_matrix_difference=comparison)
                np.savez_compressed(directory / 'cycle.npz', initial_state=x, parameters=candidate_p,
                                    period=t, weights=weights, quotient_monodromy=mq, contrast_monodromies=blocks)
                item['cycle_sha256'] = sha(directory / 'cycle.npz')
                state, period, p = x, t, candidate_p
                attempted_seed = True
                step = min(.005, step * 1.5)
            except Exception as error:
                item.update(status='failed', error=repr(error))
                if not attempted_seed or step <= .0003125:
                    record['status'] = 'stopped after bounded failure'
                else:
                    step /= 2
            item['wall_seconds'] = time.perf_counter() - started
            write(directory / 'analysis.json', item)
            record['attempts'].append(item)
            record.update(wall_seconds=time.perf_counter() - began, final_successful_rho=float(p[0]))
            write(out / 'analysis.json', record)
            print(json.dumps(dict(index=attempt, rho=rho, status=item['status'], period=item.get('period'),
                                  quotient_max=item.get('spectrum', {}).get('quotient_max_nonflow'),
                                  transverse=item.get('spectrum', {}).get('contrast_maxima'))), flush=True)
            if record['status'] != 'running':
                break
            if rho == target and item['status'] == 'complete':
                record['status'] = 'complete'
                break
        if record['status'] == 'running':
            record['status'] = 'maximum attempts incomplete'
    except Exception as error:
        record.update(status='error', error=repr(error))
        raise
    finally:
        record['wall_seconds'] = time.perf_counter() - began
        write(out / 'analysis.json', record)


if __name__ == '__main__':
    main()

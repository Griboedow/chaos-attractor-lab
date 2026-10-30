"""Reintegrate the archived original row-46 start in fixed-pair log coordinates.

This is a precision diagnostic, not a basin or Lyapunov estimate. The archived
initial state is used before any coalescence; no captured endpoint is repaired.
"""
import argparse
import datetime
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

import motor_rows as rows
import motor_pair_coordinates as pair
import motor_log_contrasts as log
import motor_log_scalar as scalar

ROOT = pair.ROOT
BASE = ROOT/'results/motor-original-seed-log'
PERM = np.array([0, 5, 1, 2, 3, 6, 4, 7])
INVERSE = np.argsort(PERM)
ORIGINAL_PAIRS = [[0, 5], [1, 2], [3, 6], [4, 7]]


def physical(z):
    x = np.concatenate([z[..., 16:23], -np.sum(z[..., 16:23], axis=-1, keepdims=True)], axis=-1)
    return np.stack([z[..., :8], z[..., 8:16], x], axis=-1)


def pack(q):
    return np.concatenate([q[..., :, 0], q[..., :, 1], q[..., :7, 2]], axis=-1)


def prepare():
    folder = BASE/'initial'
    folder.mkdir(parents=True, exist_ok=False)
    source = ROOT/'results/motor-rows/activity-path/N8.npz'
    archived = np.load(source)
    original, p = archived['initial_states'][46].copy(), archived['parameters'][46].copy()
    reordered = pack(physical(original)[PERM])
    chart = pair.to_chart(reordered)
    z = log.to_log(chart)
    assert np.all(chart[11:19] == 0.) and np.all(chart[19:23] != 0.)
    errors = dict(reordered_constraint_roundtrip=float(np.max(abs(physical(reordered)[INVERSE]-physical(original)))),
                  chart_roundtrip_original_labels=float(np.max(abs(pair.chart_rows(chart)[INVERSE]-physical(original)))),
                  log_roundtrip_original_labels=float(np.max(abs(pair.chart_rows(log.to_chart(z))[INVERSE]-physical(original)))))
    assert max(errors.values()) < 1e-15
    np.savez_compressed(folder/'trajectory.npz', parameters=p, original_initial_state=original,
                        reordered_initial_state=reordered, final_chart_state=chart,
                        final_log_state=z, time=np.array([0.]), permutation=PERM)
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='prepared',
                  source=source.relative_to(ROOT).as_posix(), source_sha256=pair.sha(source), row_index=46,
                  parameters=p.tolist(), original_initial_state=original.tolist(),
                  exact_archived_initial_state_copy=original.tobytes() == archived['initial_states'][46].tobytes(),
                  permutation_current_from_original=PERM.tolist(), original_pairs=ORIGINAL_PAIRS,
                  floating_conversion_errors=errors, initial_log10_halfcontrast_norms=log.log10_norms(z).tolist(),
                  scalar_algebra=log.algebra_audit(p, scalar.rhs), chart_algebra=pair.audit(p),
                  no_captured_endpoint_or_new_seed=True, code_sha256=pair.sha(__file__))
    pair.write(folder/'initial-audit.json', record)
    (folder/'source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(record, indent=2))


def short_control():
    folder = BASE/'initial'
    data = np.load(folder/'trajectory.npz')
    p, chart, z, original = [data[k] for k in ['parameters', 'final_chart_state', 'final_log_state', 'original_initial_state']]
    checkpoints = np.linspace(0., 20., 201)
    common = dict(method='DOP853', max_step=.025, rtol=2e-12, atol=2e-14, t_eval=checkpoints)
    cart = solve_ivp(lambda t, q: rows.rhs(t, q, p, 8), (0., 20.), original, **common)
    direct = solve_ivp(lambda t, q: pair.rhs(t, q, p), (0., 20.), chart, **common)
    logarithmic = solve_ivp(lambda t, q: scalar.rhs(t, q, p), (0., 20.), z, **common)
    assert cart.success and direct.success and logarithmic.success
    cart_rows = physical(cart.y.T)
    chart_rows = pair.chart_rows(direct.y.T)[:, INVERSE]
    log_rows = pair.chart_rows(log.to_chart(logarithmic.y.T))[:, INVERSE]
    errors = dict(cartesian_chart=float(np.max(abs(cart_rows-chart_rows))),
                  cartesian_log=float(np.max(abs(cart_rows-log_rows))),
                  chart_log=float(np.max(abs(chart_rows-log_rows))))
    record = dict(status='passed' if max(errors.values()) < 1e-8 else 'failed', duration=20.,
                  maximum_row_triple_errors=errors, threshold=1e-8,
                  source_data_sha256=pair.sha(folder/'trajectory.npz'),
                  nfev=[cart.nfev, direct.nfev, logarithmic.nfev],
                  scope='Original archived Cartesian start compared with converted chart/log starts in original row labels.')
    pair.write(folder/'short-control.json', record)
    if record['status'] != 'passed':
        raise ArithmeticError(record)
    print(json.dumps(record))


def run(args):
    folder = BASE/args.name
    folder.mkdir(parents=True, exist_ok=False)
    source = BASE/'initial/trajectory.npz'
    data = np.load(source)
    p, initial = data['parameters'].copy(), data['final_log_state'].copy()
    audit_file = BASE/'scalar-audit-20tau/audit.json'
    audit = json.loads(audit_file.read_text())
    assert audit['status'] == 'passed' and audit['code_sha256'] == pair.sha(scalar.__file__)
    assert audit['log_code_sha256'] == pair.sha(log.__file__)
    assert audit['source_data_sha256'] == pair.sha(source)
    assert json.loads((BASE/'initial/short-control.json').read_text())['status'] == 'passed'
    tracker, times, states = log.UnderflowTracker(), [0.], [initial.copy()]
    summary = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='running',
                   parameters=p.tolist(), horizon=args.horizon, max_step=args.max_step,
                   rtol=args.rtol, atol=args.atol, original_pairs=ORIGINAL_PAIRS,
                   checkpoint_spacing=1., source_data_sha256=pair.sha(source),
                   scalar_audit_sha256=pair.sha(audit_file), code_sha256=pair.sha(__file__),
                   no_noise_projection_or_renormalization=True, blocks=[],
                   limitation='Floating-point means and directions; independent runs may separate in sensitive dynamics.')
    began = time.perf_counter()
    for left in np.arange(0., args.horizon, 500.):
        right = min(float(left+500.), args.horizon)
        checkpoints = np.arange(left+1., right+1.)
        started = time.perf_counter()
        sol = solve_ivp(lambda t, q: scalar.rhs(t, q, p, tracker), (left, right), states[-1],
                        method='DOP853', max_step=args.max_step, rtol=args.rtol,
                        atol=np.full(27, args.atol), t_eval=checkpoints)
        if not sol.success:
            summary['status'], summary['failure'] = 'failed', sol.message
            pair.write(folder/'summary.json', summary)
            raise RuntimeError(sol.message)
        times.extend(sol.t.tolist())
        states.extend(sol.y.T.copy())
        values = np.asarray(states)
        assert np.all(np.isfinite(values)) and np.all(log.gauge_norm_squared(values) > 0)
        norms = log.log10_norms(values)
        summary['blocks'].append(dict(end_tau=right, wall_seconds=time.perf_counter()-started, nfev=sol.nfev,
                                      endpoint_log10_norms=norms[-1].tolist()))
        summary.update(elapsed_tau=right, wall_seconds=time.perf_counter()-began,
                       minimum_log10_norms=np.min(norms, axis=0).tolist(), final_log10_norms=norms[-1].tolist(),
                       maximum_gauge_norm_squared_drift=float(np.max(abs(log.gauge_norm_squared(values)-log.gauge_norm_squared(initial)))),
                       underflow=tracker.records, status='complete' if right == args.horizon else 'running')
        np.savez_compressed(folder/'trajectory.npz', time=np.asarray(times), log_states=values,
                            log10_halfcontrast_norms=norms, parameters=p, initial_log_state=initial,
                            final_log_state=values[-1], permutation=PERM, original_initial_state=data['original_initial_state'])
        summary['trajectory_sha256'] = pair.sha(folder/'trajectory.npz')
        pair.write(folder/'summary.json', summary)
        print(json.dumps(summary['blocks'][-1]), flush=True)
    for name, file in [('source.py', __file__), ('scalar-source.py', scalar.__file__),
                       ('log-source.py', log.__file__), ('pair-source.py', pair.__file__)]:
        (folder/name).write_bytes(Path(file).read_bytes())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['prepare', 'short-control', 'run'])
    parser.add_argument('--output', type=Path, default=BASE,
                        help='Experiment root; use a fresh outputs/ subfolder for reproduction.')
    parser.add_argument('--name', default='step-0p1')
    parser.add_argument('--horizon', type=float, default=10000.)
    parser.add_argument('--max-step', type=float, default=.1)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--atol', type=float, default=1e-12)
    args = parser.parse_args()
    BASE = args.output.resolve()
    {'prepare': prepare, 'short-control': short_control, 'run': lambda: run(args)}[args.mode]()

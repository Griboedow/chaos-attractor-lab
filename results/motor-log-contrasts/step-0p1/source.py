"""Logarithmic gauge coordinates for nonzero fixed-pair motor contrasts.

The exact 27-state equations represent the 23 physical states plus four constant
gauge norms. Directions are never renormalized and contrasts are never projected.
Floating-point underflow of quadratic mean feedback is recorded explicitly.
"""

import argparse
import datetime
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import linear_sum_assignment, minimize_scalar

import motor_pair_coordinates as pair

ROOT = pair.ROOT
LN10 = math.log(10.)


def to_log(chart):
    chart = np.asarray(chart)
    norms = pair.contrast_norms(chart)
    if np.any(norms <= 0):
        raise ValueError('Log coordinates require four nonzero pair contrasts')
    return np.concatenate([chart[..., :11], chart[..., 11:15]/norms,
                           chart[..., 15:19]/norms, chart[..., 19:23]/norms,
                           np.log(norms)], axis=-1)


def to_chart(z):
    """Lossy when exp(ell) underflows; never used to advance the log system."""
    factor = np.exp(z[..., 23:27])
    return np.concatenate([z[..., :11], z[..., 11:15]*factor,
                           z[..., 15:19]*factor, z[..., 19:23]*factor], axis=-1)


def gauge_norm_squared(z):
    return z[..., 11:15]**2+z[..., 15:19]**2+z[..., 19:23]**2


def log10_norms(z):
    return (z[..., 23:27]+.5*np.log(gauge_norm_squared(z)))/LN10


def rhs(t, z, p, tracker=None):
    rho, r, b = p
    c, s = z[:4], z[4:8]
    x = np.r_[z[8:11], -np.sum(z[8:11])]
    uc, us, ux, ell = z[11:15], z[15:19], z[19:23], z[23:27]
    v = -rho*x+r*(s-np.mean(s))
    duv = -rho*ux+r*us
    auc = -uc-v*us-s*duv
    aus = -us+v*uc+c*duv
    aux = duv
    normsq = uc*uc+us*us+ux*ux
    alpha = (uc*auc+us*aus+ux*aux)/normsq
    weight = np.exp(2*ell)
    f = np.empty_like(z)
    f[:4] = 1-c-v*s-weight*duv*us
    f[4:8] = b-s+v*c+weight*duv*uc
    f[8:11] = v[:3]
    f[11:15] = auc-alpha*uc
    f[15:19] = aus-alpha*us
    f[19:23] = aux-alpha*ux
    f[23:27] = alpha
    if tracker is not None:
        tracker.observe(t, z, weight, p)
    return f


class UnderflowTracker:
    def __init__(self):
        self.records = [dict(pair=label, first_rhs_evaluation_tau=None,
                             number_of_zero_weight_rhs_evaluations=0,
                             maximum_omitted_component_log_bound=None,
                             maximum_two_ell_when_weight_zero=None) for label in pair.PAIRS]

    def observe(self, t, z, weight, p):
        for j in np.flatnonzero(weight == 0):
            record = self.records[j]
            record['number_of_zero_weight_rhs_evaluations'] += 1
            if record['first_rhs_evaluation_tau'] is None:
                record['first_rhs_evaluation_tau'] = float(t)
            umax = max(abs(z[11+j]), abs(z[15+j]), abs(z[19+j]))
            bound = float(2*z[23+j]+math.log(abs(p[0])+abs(p[1]))+2*math.log(umax))
            previous = record['maximum_omitted_component_log_bound']
            record['maximum_omitted_component_log_bound'] = bound if previous is None else max(previous, bound)
            previous = record['maximum_two_ell_when_weight_zero']
            record['maximum_two_ell_when_weight_zero'] = float(2*z[23+j]) if previous is None else max(previous, float(2*z[23+j]))


def algebra_audit(p):
    rng = np.random.default_rng(2026103052)
    errors = dict(field_scaled=0., gauge_derivative=0., reconstruction_scaled=0., gauge_equivalence_scaled=0.)
    for _ in range(32):
        chart = rng.normal(size=23)
        z = to_log(chart)
        # Arbitrary nonunit gauge norms must obey the same physical equations.
        gauge = np.exp(rng.uniform(-2, 2, 4))
        z[11:15] *= gauge
        z[15:19] *= gauge
        z[19:23] *= gauge
        z[23:27] -= np.log(gauge)
        f = rhs(0., z, p)
        factor = np.exp(z[23:27])
        pushed = np.r_[f[:11], factor*(f[11:15]+z[11:15]*f[23:27]),
                       factor*(f[15:19]+z[15:19]*f[23:27]),
                       factor*(f[19:23]+z[19:23]*f[23:27])]
        expected = pair.rhs(0., chart, p)
        errors['field_scaled'] = max(errors['field_scaled'], float(np.max(abs(pushed-expected))/(1+np.max(abs(expected)))))
        norm_derivative = 2*(z[11:15]*f[11:15]+z[15:19]*f[15:19]+z[19:23]*f[19:23])
        errors['gauge_derivative'] = max(errors['gauge_derivative'], float(np.max(abs(norm_derivative))/(1+np.max(gauge_norm_squared(z)))))
        errors['reconstruction_scaled'] = max(errors['reconstruction_scaled'], float(np.max(abs(to_chart(z)-chart))/(1+np.max(abs(chart)))))
        plain = rhs(0., to_log(chart), p)
        errors['gauge_equivalence_scaled'] = max(errors['gauge_equivalence_scaled'], float(np.max(abs(plain[:11]-f[:11]))/(1+np.max(abs(f[:11])))))
    if max(errors.values()) > 1e-12:
        raise ArithmeticError(errors)
    return dict(samples=32, arbitrary_nonunit_gauges=True, errors=errors,
                exact_identity='D(reconstruction) F_log = F_chart; d||u_a||^2/dtau=0',
                no_renormalization_or_projection=True)


def short_control(chart, p, start, duration=200.):
    times = np.linspace(start, start+duration, 201)
    common = dict(method='DOP853', rtol=2e-12, max_step=.025, t_eval=times)
    a = solve_ivp(lambda t, q: pair.rhs(t, q, p), (start, start+duration), chart,
                  atol=np.r_[np.full(11, 2e-14), np.full(12, 1e-200)], **common)
    b = solve_ivp(lambda t, q: rhs(t, q, p), (start, start+duration), to_log(chart),
                  atol=np.full(27, 2e-14), **common)
    if not(a.success and b.success):
        raise RuntimeError([a.message, b.message])
    ac, bl = a.y.T, b.y.T
    differences = dict(maximum_mean_difference=float(np.max(abs(ac[:, :11]-bl[:, :11]))),
                       maximum_chart_difference=float(np.max(abs(ac-to_chart(bl)))),
                       maximum_log10_norm_difference=float(np.max(abs(np.log10(pair.contrast_norms(ac))-log10_norms(bl)))),
                       maximum_gauge_norm_squared_drift=float(np.max(abs(gauge_norm_squared(bl)-gauge_norm_squared(bl[0])))))
    if max(differences.values()) > 1e-7:
        raise ArithmeticError(differences)
    return dict(start_tau=start, duration=duration, checkpoints=len(times), differences=differences,
                acceptance_threshold=1e-7, chart_nfev=a.nfev, log_nfev=b.nfev,
                scope='Short comparison at representable contrasts; not a long-time error bound.')


def make_events():
    events, descriptions = [], []
    for j in range(4):
        def event(t, z, j=j):
            return z[23+j]+.5*math.log(z[11+j]**2+z[15+j]**2+z[19+j]**2)+4*LN10
        event.direction, event.terminal = 1, False
        events.append(event)
        descriptions.append(dict(pair_index=j, pair=pair.PAIRS[j], threshold_log10=-4, times=[]))
    return events, descriptions


def run(args):
    folder, parent = args.output.resolve(), args.resume_from.resolve()
    if folder.exists():
        raise FileExistsError(folder)
    saved = np.load(parent/'trajectory.npz')
    meta = json.loads((parent/'summary.json').read_text())
    assert meta['status'] == 'complete'
    p, start = saved['parameters'], float(saved['time'][-1])
    chart = saved['final_chart_state'].copy()
    assert np.array_equal(chart, saved['chart_states'][-1])
    initial = z = to_log(chart)
    if args.horizon <= start:
        raise ValueError('Horizon must exceed the saved chart time')
    folder.mkdir(parents=True)
    (folder/'source.py').write_bytes(Path(__file__).read_bytes())
    (folder/'pair-source.py').write_bytes(Path(pair.__file__).read_bytes())
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='initializing',
                  parent=parent.relative_to(ROOT).as_posix(), parent_data_sha256=pair.sha(parent/'trajectory.npz'),
                  parent_summary_sha256=pair.sha(parent/'summary.json'), code_sha256=pair.sha(__file__),
                  pair_code_sha256=pair.sha(pair.__file__), parameters=p.tolist(), start_tau=start,
                  reference_period=meta['reference_period'], horizon_requested=args.horizon,
                  max_step=args.max_step, rtol=args.rtol, mean_atol=args.mean_atol,
                  direction_atol=args.direction_atol, logarithm_atol=args.logarithm_atol,
                  coordinate_order='mean_c4,mean_s4,mean_X3,u_c4,u_s4,u_X4,ell4',
                  coordinate_map='physical halfcontrast delta_a=exp(ell_a)*u_a',
                  no_noise_projection_or_renormalization=True,
                  initial_chart_reconstruction_max_error=float(np.max(abs(to_chart(initial)-chart))),
                  gauge_norm_squared_initial=gauge_norm_squared(initial).tolist(),
                  underflow_scope='Zero exp(2ell) RHS stage evaluations, including trial/rejected stages. Bounds concern instantaneous omitted mean forcing, not accumulated trajectory error.',
                  blocks=[], limitations=[
                      'The exact algebra is a coordinate description; floating-point integration is not a proof.',
                      'Quadratic feedback can underflow or be rounded away before affecting O(1) means.',
                      'Four auxiliary gauge invariants are not physical neutral modes.',
                      'No heteroclinic-connection, asymptotic attraction, or chaos theorem follows.'])
    pair.write(folder/'summary.json', record)
    record['algebra'] = algebra_audit(p)
    record['short_equivalence_control'] = short_control(chart, p, start)
    events, descriptions = make_events()
    record['events'] = descriptions
    tracker = UnderflowTracker()
    record['status'] = 'running'
    period, current = meta['reference_period'], start
    times, states = [start], [initial.copy()]
    began = time.perf_counter()
    atol = np.r_[np.full(11, args.mean_atol), np.full(12, args.direction_atol), np.full(4, args.logarithm_atol)]
    while current < args.horizon:
        block_start = current
        end = min(args.horizon, current+args.chunk)
        first = int(math.floor(current/period))+1
        checkpoints = np.arange(first, int(math.floor(end/period))+1)*period
        checkpoints = np.r_[checkpoints[(checkpoints > current) & (checkpoints < end)], end]
        began_block = time.perf_counter()
        sol = solve_ivp(lambda t, q: rhs(t, q, p, tracker), (current, end), z,
                        method='DOP853', max_step=args.max_step, rtol=args.rtol, atol=atol,
                        t_eval=checkpoints, events=events)
        if not sol.success or not np.all(np.isfinite(sol.y)):
            record.update(status='failed', error=sol.message)
            break
        block = sol.y.T
        if np.any(gauge_norm_squared(block) <= 0):
            record.update(status='failed', error='Zero gauge direction')
            break
        for j, event_times in enumerate(sol.t_events):
            descriptions[j]['times'].extend(event_times.tolist())
        times.extend(sol.t.tolist())
        states.extend(block.copy())
        z, current = block[-1].copy(), end
        logs = log10_norms(block)
        record['blocks'].append(dict(start_tau=block_start,
                                     first_checkpoint_tau=float(sol.t[0]), end_tau=end, nfev=sol.nfev,
                                     wall_seconds=time.perf_counter()-began_block,
                                     minimum_log10_halfcontrast_norms=logs.min(axis=0).tolist(),
                                     endpoint_log10_halfcontrast_norms=logs[-1].tolist(),
                                     maximum_gauge_norm_squared_drift=float(np.max(abs(gauge_norm_squared(block)-gauge_norm_squared(initial))))))
        stored = np.array(states)
        record.update(elapsed_tau=current, wall_seconds=time.perf_counter()-began,
                      underflow=tracker.records, final_log_state=z.tolist(),
                      maximum_gauge_norm_squared_drift=float(np.max(abs(gauge_norm_squared(stored)-gauge_norm_squared(initial)))))
        np.savez_compressed(folder/'trajectory.npz', parameters=p, time=np.array(times), log_states=stored,
                            log10_halfcontrast_norms=log10_norms(stored), initial_chart_state=chart,
                            initial_log_state=initial, final_log_state=z,
                            origin_initial_chart_state=saved['origin_initial_chart_state'])
        pair.write(folder/'summary.json', record)
        print(json.dumps(dict(tau=end, wall=round(record['wall_seconds'], 2),
                              log10_halfcontrast=logs[-1].tolist(),
                              underflow_counts=[r['number_of_zero_weight_rhs_evaluations'] for r in tracker.records])), flush=True)
    if record['status'] == 'running':
        record['status'] = 'complete'
    if (folder/'trajectory.npz').exists():
        record['trajectory_sha256'] = pair.sha(folder/'trajectory.npz')
    pair.write(folder/'summary.json', record)


def analyze(args):
    folder = args.output.resolve()
    if (folder/'diagnostics.json').exists():
        raise FileExistsError(folder/'diagnostics.json')
    data = np.load(folder/'trajectory.npz')
    meta = json.loads((folder/'summary.json').read_text())
    times, states, logs = data['time'], data['log_states'], data['log10_halfcontrast_norms']
    recovered, pending = [], []
    for j in range(4):
        last = -1.
        for index in np.flatnonzero(logs[:, j] < -12):
            if times[index] <= last:
                continue
            future = [t for t in meta['events'][j]['times'] if t > times[index]]
            if not future:
                pending.append(dict(pair=pair.PAIRS[j], first_low_checkpoint_tau=float(times[index]),
                                    minimum_log10_halfcontrast=float(np.min(logs[index:, j]))))
                break
            departure = future[0]
            mask = (times >= times[index]) & (times <= departure)
            recovered.append(dict(pair=pair.PAIRS[j], first_low_checkpoint_tau=float(times[index]),
                                  minimum_log10_halfcontrast=float(np.min(logs[mask, j])),
                                  recovery_upward_crossing_1e_minus_4=departure))
            last = departure
    recovered.sort(key=lambda r: r['recovery_upward_crossing_1e_minus_4'])
    ref = np.load(args.source)
    reference = pair.physical(ref['full_states'])
    reference[-1] = reference[0]
    period = float(ref['time'][-1])
    spline = CubicSpline(ref['time'], reference, bc_type='periodic')
    coarse = np.arange(80)*period/80
    coarse_rows = spline(coarse)
    matches = []
    for index in range(0, len(times), 10):
        current = pair.chart_rows(to_chart(states[index]))
        def cost(phase, details=False):
            rows = spline(phase % period)
            sq = np.sum((current[:, None, :]-rows[None, :, :])**2, axis=2)
            ri, ci = linear_sum_assignment(sq)
            value = float(sq[ri, ci].sum())
            return (math.sqrt(value), ci) if details else value
        values = []
        for rows in coarse_rows:
            sq = np.sum((current[:, None, :]-rows[None, :, :])**2, axis=2)
            ri, ci = linear_sum_assignment(sq)
            values.append(sq[ri, ci].sum())
        guess = coarse[int(np.argmin(values))]
        fit = minimize_scalar(cost, bounds=(guess-period/80, guess+period/80), method='bounded', options={'xatol': 1e-12})
        distance, perm = cost(fit.x, True)
        matches.append(dict(tau=float(times[index]), distance=distance, phase=float(fit.x % period),
                            permutation_current_from_reference=perm.tolist(),
                            contracting_pair=np.flatnonzero(np.isin(perm, [2, 3])).tolist(),
                            expanding_pair=np.flatnonzero(np.isin(perm, [6, 7])).tolist(),
                            halfcontrast_log10=logs[index].tolist()))
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  minimum_log10_halfcontrast_norms=logs.min(axis=0).tolist(),
                  final_log10_halfcontrast_norms=logs[-1].tolist(), recovered_tiny_contrasts=recovered,
                  unresolved_small_contrasts_at_end=pending, orbit_matches=matches,
                  recovery_criterion='After checkpoint below1e-12, first subsequent upward crossing of1e-4',
                  finite_log_states=bool(np.all(np.isfinite(states))),
                  positive_gauge_norms=bool(np.all(gauge_norm_squared(states) > 0)),
                  data_sha256=pair.sha(folder/'trajectory.npz'), code_sha256=pair.sha(__file__),
                  caveat='Reconstructed rows used for orbit matching discard tiny contrasts; matching is a coarse dynamical diagnostic with interpolation and floating-point floors.')
    pair.write(folder/'diagnostics.json', record)
    (folder/'analysis.source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps({k: record[k] for k in ['minimum_log10_halfcontrast_norms', 'final_log10_halfcontrast_norms', 'recovered_tiny_contrasts']}))


def compare(args):
    folders = [args.output.resolve(), args.other.resolve()]
    data = [np.load(f/'trajectory.npz') for f in folders]
    a, b = data
    assert np.array_equal(a['time'], b['time'])
    assert np.array_equal(a['origin_initial_chart_state'], b['origin_initial_chart_state'])
    assert np.array_equal(a['parameters'], b['parameters'])
    events = []
    diagnostics = [json.loads((f/'diagnostics.json').read_text()) for f in folders]
    for event in diagnostics[0]['recovered_tiny_contrasts']:
        candidates = [q for q in diagnostics[1]['recovered_tiny_contrasts'] if q['pair'] == event['pair']]
        if candidates:
            other = min(candidates, key=lambda q: abs(q['recovery_upward_crossing_1e_minus_4']-event['recovery_upward_crossing_1e_minus_4']))
            ta, tb = event['recovery_upward_crossing_1e_minus_4'], other['recovery_upward_crossing_1e_minus_4']
            events.append(dict(pair=event['pair'], first_tau=ta, second_tau=tb, difference=tb-ta))
    windows = []
    for end in [50000., 75000., 100000., 125000., 150000.]:
        mask = a['time'] <= end
        if not np.any(mask) or end > a['time'][-1]:
            continue
        windows.append(dict(through_tau=end,
                            maximum_mean_difference=float(np.max(abs(a['log_states'][mask, :11]-b['log_states'][mask, :11]))),
                            maximum_log10_norm_difference=float(np.max(abs(a['log10_halfcontrast_norms'][mask]-b['log10_halfcontrast_norms'][mask])))))
    record = dict(folders=[f.relative_to(ROOT).as_posix() for f in folders],
                  exactly_equal_origin_chart_states=True, exactly_equal_checkpoint_times=True,
                  initial_log_state_difference=float(np.max(abs(a['initial_log_state']-b['initial_log_state']))),
                  windows=windows, recovery_event_comparison=events,
                  data_sha256=[pair.sha(f/'trajectory.npz') for f in folders], code_sha256=pair.sha(__file__),
                  caveat='Independent histories include the preceding23-state controls. Agreement is empirical, not an error bound for the exact far-time orbit.')
    pair.write(folders[0]/'comparison.json', record)
    (folders[0]/'comparison.source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['run', 'analyze', 'compare'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume-from', type=Path)
    parser.add_argument('--other', type=Path)
    parser.add_argument('--source', type=Path, default=ROOT/'results/motor-postcritical/six-cluster-cycle-rho-3p268/cycle.npz')
    parser.add_argument('--horizon', type=float, default=150000.)
    parser.add_argument('--chunk', type=float, default=500.)
    parser.add_argument('--max-step', type=float, default=.1)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--mean-atol', type=float, default=1e-12)
    parser.add_argument('--direction-atol', type=float, default=1e-12)
    parser.add_argument('--logarithm-atol', type=float, default=1e-12)
    args = parser.parse_args()
    {'run': run, 'analyze': analyze, 'compare': compare}[args.command](args)

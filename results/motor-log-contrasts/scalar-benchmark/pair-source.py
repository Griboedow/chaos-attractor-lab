"""Exact mean/half-difference chart for the full eight-row motor equations.

This precision diagnostic keeps every one of the 23 states. No noise, cluster
projection or changed feedback is applied during integration.
"""

import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import linear_sum_assignment, minimize_scalar

from motor_rows import field_jac as full_field_jac, rhs as full_rhs

ROOT = Path(__file__).resolve().parents[1]
PAIRS = [[0, 1], [2, 3], [4, 5], [6, 7]]


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def chart_rows(z):
    z = np.asarray(z)
    means = np.stack([z[..., :4], z[..., 4:8],
                      np.concatenate([z[..., 8:11], -np.sum(z[..., 8:11], axis=-1, keepdims=True)], axis=-1)], axis=-1)
    diffs = np.stack([z[..., 11:15], z[..., 15:19], z[..., 19:23]], axis=-1)
    return np.stack([means + diffs, means - diffs], axis=-2).reshape(z.shape[:-1] + (8, 3))


def to_full(z):
    q = chart_rows(z)
    return np.concatenate([q[..., :, 0], q[..., :, 1], q[..., :7, 2]], axis=-1)


def to_chart(z):
    z = np.asarray(z)
    x = np.concatenate([z[..., 16:23], -np.sum(z[..., 16:23], axis=-1, keepdims=True)], axis=-1)
    q = np.stack([z[..., :8], z[..., 8:16], x], axis=-1)
    means = (q[..., ::2, :] + q[..., 1::2, :]) / 2
    diffs = (q[..., ::2, :] - q[..., 1::2, :]) / 2
    return np.concatenate([means[..., :, 0], means[..., :, 1], means[..., :3, 2],
                           diffs[..., :, 0], diffs[..., :, 1], diffs[..., :, 2]], axis=-1)


E = np.column_stack([to_full(v) for v in np.eye(23)])
EINV = np.column_stack([to_chart(v) for v in np.eye(23)])


def rhs(t, z, p):
    rho, r, b = p
    c, s = z[:4], z[4:8]
    x = np.empty(4, dtype=z.dtype)
    x[:3], x[3] = z[8:11], -np.sum(z[8:11])
    dc, ds, dx = z[11:15], z[15:19], z[19:23]
    v = -rho*x + r*(s-np.mean(s))
    dv = -rho*dx + r*ds
    f = np.empty_like(z)
    f[:4] = 1-c-v*s-dv*ds
    f[4:8] = b-s+v*c+dv*dc
    f[8:11] = v[:3]
    f[11:15] = -dc-v*ds-s*dv
    f[15:19] = -ds+v*dc+c*dv
    f[19:23] = dv
    return f


def contrast_norms(z):
    """Scaled norms avoid squaring tiny contrasts into floating-point zero."""
    d = np.stack([z[..., 11:15], z[..., 15:19], z[..., 19:23]], axis=-1)
    scale = np.max(abs(d), axis=-1)
    safe = np.where(scale == 0, 1., scale)
    return scale * np.sqrt(np.sum((d/safe[..., None])**2, axis=-1))


def one_norm(z, pair):
    a, b, c = z[11+pair], z[15+pair], z[19+pair]
    scale = max(abs(a), abs(b), abs(c))
    return scale * math.sqrt((a/scale)**2 + (b/scale)**2 + (c/scale)**2) if scale else 0.


def audit(p):
    assert np.array_equal(EINV @ E, np.eye(23))
    assert np.array_equal(E @ EINV, np.eye(23))
    rng = np.random.default_rng(2026103051)
    errors = dict(field_scaled=0., jacobian_scaled=0., inverse_roundtrip=0.)
    for _ in range(32):
        z = rng.normal(size=23)
        full, jf = full_field_jac(to_full(z), p, 8)
        f = rhs(0., z, p)
        jc = np.column_stack([np.imag(rhs(0., z.astype(complex)+1e-28j*v, p))/1e-28 for v in np.eye(23)])
        errors['field_scaled'] = max(errors['field_scaled'], float(np.max(abs(E@f-full))/(1+np.max(abs(full)))))
        expected = EINV@jf@E
        errors['jacobian_scaled'] = max(errors['jacobian_scaled'], float(np.max(abs(jc-expected))/(1+np.max(abs(expected)))))
        errors['inverse_roundtrip'] = max(errors['inverse_roundtrip'], float(np.max(abs(to_chart(to_full(z))-z))))
    if max(errors.values()) > 1e-12:
        raise ArithmeticError(errors)
    z = rng.normal(size=23)
    zero_checks = []
    for pair in range(4):
        indices = [11+pair, 15+pair, 19+pair]
        z[indices] = 0.
        zero_checks.append(bool(np.all(rhs(0., z, p)[indices] == 0)))
    assert all(zero_checks)
    z[11:15], z[15:19], z[19:23] = 1e-80, -2e-80, 3e-80
    tiny = contrast_norms(z)
    rebuilt = contrast_norms(to_chart(to_full(z)))
    return dict(samples=32, errors=errors, exact_matrix_inverses=True,
                homogeneous_zero_contrasts=zero_checks,
                tiny_chart_norms=tiny.tolist(), tiny_reconstructed_cartesian_norms=rebuilt.tolist(),
                meaning='The rebuilt Cartesian rows lose tiny contrasts; the direct chart retains them.')


def initialize(source, unstable=1e-7, seed=1e-12, sign=1):
    a = np.load(source)
    p, period = a['parameters'], float(a['time'][-1])
    z = to_chart(a['full_states'][0])
    monodromy = EINV @ a['full_monodromy'] @ E
    modes = []
    for pair, amplitude in [(3, sign*unstable), (1, seed)]:
        indices = [11+pair, 15+pair, 19+pair]
        block = monodromy[np.ix_(indices, indices)]
        eig, vectors = np.linalg.eig(block)
        index = int(np.argmax(abs(eig)))
        assert abs(eig[index].imag) < 1e-10
        vector = vectors[:, index].real
        vector /= np.linalg.norm(vector)
        if vector[0] < 0:
            vector *= -1
        before = z[indices].copy()
        z[indices] = amplitude*vector
        leakage = monodromy[:, indices].copy()
        leakage[indices, :] = 0.
        modes.append(dict(pair=PAIRS[pair], pair_index=pair, halfcontrast_amplitude=amplitude,
                          block=block.tolist(), eigenvalues=[dict(real=float(v.real), imag=float(v.imag)) for v in eig],
                          selected_multiplier=float(eig[index].real),
                          exponent=float(np.log(abs(eig[index]))/period), direction=vector.tolist(),
                          archived_halfcontrast=before.tolist(), off_block_max=float(np.max(abs(leakage)))))
    return z, p, period, modes


def short_control(z, p, period):
    kwargs = dict(method='DOP853', rtol=2e-12, max_step=.025)
    atol = np.r_[np.full(11, 2e-14), np.full(12, 1e-150)]
    chart = solve_ivp(lambda t, q: rhs(t, q, p), (0., 2*period), z, atol=atol, **kwargs)
    full = solve_ivp(lambda t, q: full_rhs(t, q, p, 8), (0., 2*period), to_full(z), atol=2e-14, **kwargs)
    if not (chart.success and full.success):
        raise RuntimeError((chart.message, full.message))
    error = float(np.max(abs(to_full(chart.y[:, -1])-full.y[:, -1])))
    if error > 1e-8:
        raise ArithmeticError(f'Short Cartesian control discrepancy {error} exceeds1e-8')
    return dict(horizon=2*period, maximum_endpoint_state_error=error, acceptance_threshold=1e-8,
                chart_nfev=chart.nfev, full_nfev=full.nfev,
                scope='Short full Cartesian control; not a long-time tiny-contrast accuracy test.')


def make_events():
    functions, descriptions = [], []
    for pair in range(4):
        for threshold in [1e-4, .01]:
            def event(t, z, pair=pair, threshold=threshold):
                return one_norm(z, pair)-threshold
            event.direction = 1
            event.terminal = False
            functions.append(event)
            descriptions.append(dict(pair=PAIRS[pair], pair_index=pair, halfcontrast_threshold=threshold,
                                     direction='upward', times=[]))
    return functions, descriptions


def run(args):
    source, output = args.source.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    z, p, period, modes = initialize(source, args.unstable, args.seed, args.sign)
    start_tau, origin_initial, resume_record = 0., z.copy(), None
    if args.resume_from:
        parent = args.resume_from.resolve()
        saved = np.load(parent/'trajectory.npz')
        previous = json.loads((parent/'summary.json').read_text())
        assert previous['status'] == 'complete'
        assert np.array_equal(saved['parameters'], p)
        assert previous['source_sha256'] == sha(source)
        z = saved['final_chart_state'].copy()
        assert np.array_equal(z, saved['chart_states'][-1])
        start_tau = float(saved['time'][-1])
        origin_initial = saved['origin_initial_chart_state'].copy() if 'origin_initial_chart_state' in saved else saved['initial_chart_state'].copy()
        resume_record = dict(parent=parent.relative_to(ROOT).as_posix(),
                             trajectory_sha256=sha(parent/'trajectory.npz'), parent_summary_sha256=sha(parent/'summary.json'),
                             start_tau=start_tau, exact_saved_state_restored=True,
                             parent_code_sha256=previous['code_sha256'])
        modes = previous['initial_modes']
        if args.horizon <= start_tau:
            raise ValueError('Continuation horizon must exceed saved final time')
    initial = z.copy()
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='initializing',
                  source=source.relative_to(ROOT).as_posix(), source_sha256=sha(source), code_sha256=sha(__file__),
                  parameters=p.tolist(), reference_period=period, pair_labels=PAIRS,
                  coordinate_order='mean_c4,mean_s4,mean_X3,half_dc4,half_ds4,half_dX4',
                  initial_modes=modes, initial_chart_state=z.tolist(),
                  horizon_requested=args.horizon, chunk_tau=args.chunk, max_step=args.max_step,
                  rtol=args.rtol, mean_atol=args.mean_atol, contrast_atol=args.contrast_atol,
                  event_definition='Upward crossings of robust Euclidean halfcontrast norms',
                  no_noise_or_projection=True, blocks=[], limitations=[
                      'Approximate archived saddle and floating-point eigendirections.',
                      'Fixed-pair chart only protects differences inside its four selected pairs.',
                      'Long-time sensitive trajectories require event and step-size comparisons.',
                      'No proof of a heteroclinic connection, chaos, or an attractor.'])
    record['resume'] = resume_record
    record['origin_initial_chart_state'] = origin_initial.tolist()
    write(output/'summary.json', record)
    (output/'source.py').write_bytes(Path(__file__).read_bytes())
    record['algebra'] = audit(p)
    record['short_full_coordinate_control'] = short_control(z, p, period) if not resume_record else dict(
        inherited_from=resume_record['parent'], reason='Continuation restores every chart coordinate; the Cartesian initial-state control is in the parent run.')
    events, descriptions = make_events()
    record['events'] = descriptions
    record['status'] = 'running'
    times, states = [start_tau], [z.copy()]
    global_t, began = start_tau, time.perf_counter()
    atol = np.r_[np.full(11, args.mean_atol), np.full(12, args.contrast_atol)]
    while global_t < args.horizon:
        end = min(args.horizon, global_t+args.chunk)
        first = int(math.floor(global_t/period))+1
        sample = np.arange(first, int(math.floor(end/period))+1)*period
        sample = sample[(sample > global_t) & (sample < end)]
        sample = np.r_[sample, end]
        wall_start = time.perf_counter()
        sol = solve_ivp(lambda t, q: rhs(t, q, p), (global_t, end), z,
                        method='DOP853', rtol=args.rtol, atol=atol, max_step=args.max_step,
                        t_eval=sample, events=events)
        if not sol.success:
            record.update(status='failed', error=sol.message)
            break
        block = sol.y.T
        norms = contrast_norms(block)
        zero_rows = np.argwhere(norms == 0)
        for j, event_times in enumerate(sol.t_events):
            descriptions[j]['times'].extend(event_times.tolist())
        times.extend(sol.t.tolist())
        states.extend(block.copy())
        z, global_t = block[-1].copy(), end
        record['blocks'].append(dict(first_checkpoint_tau=float(sol.t[0]) if len(sol.t) else global_t,
                                     end_tau=end, nfev=sol.nfev, wall_seconds=time.perf_counter()-wall_start,
                                     minimum_halfcontrast_norms=norms.min(axis=0).tolist(),
                                     maximum_halfcontrast_norms=norms.max(axis=0).tolist(),
                                     endpoint_log10_halfcontrast_norms=np.log10(np.maximum(norms[-1], np.nextafter(0., 1.))).tolist()))
        record.update(elapsed_tau=global_t, wall_seconds=time.perf_counter()-began,
                      final_chart_state=z.tolist())
        stored = np.array(states)
        stored_norms = contrast_norms(stored)
        np.savez_compressed(output/'trajectory.npz', parameters=p, time=np.array(times), chart_states=stored,
                            halfcontrast_norms=stored_norms, full_states=to_full(stored),
                            initial_chart_state=initial, final_chart_state=z, origin_initial_chart_state=origin_initial,
                            E=E, E_inverse=EINV)
        if len(zero_rows):
            record.update(status='underflow_or_exact_zero', first_zero_checkpoint_tau=float(sol.t[zero_rows[0, 0]]),
                          first_zero_pair_index=int(zero_rows[0, 1]))
        write(output/'summary.json', record)
        print(json.dumps(dict(tau=end, wall=round(record['wall_seconds'], 2),
                              log10_halfcontrast=record['blocks'][-1]['endpoint_log10_halfcontrast_norms'],
                              status=record['status'])), flush=True)
        if record['status'] != 'running':
            break
    if record['status'] == 'running':
        record['status'] = 'complete'
    record['trajectory_sha256'] = sha(output/'trajectory.npz')
    write(output/'summary.json', record)


def physical(full):
    full = np.asarray(full)
    x = np.concatenate([full[..., 16:], -np.sum(full[..., 16:], axis=-1, keepdims=True)], axis=-1)
    return np.stack([full[..., :8], full[..., 8:16], x], axis=-1)


def analyze(args):
    folder, source = args.output.resolve(), args.source.resolve()
    target = folder/'diagnostics.json'
    if target.exists():
        raise FileExistsError(target)
    data, ref = np.load(folder/'trajectory.npz'), np.load(source)
    times, states, norms = data['time'], data['chart_states'], data['halfcontrast_norms']
    rows = chart_rows(states)
    reference = physical(ref['full_states'])
    reference[-1] = reference[0]
    period = float(ref['time'][-1])
    spline = CubicSpline(ref['time'], reference, bc_type='periodic')
    coarse = np.arange(80)*period/80
    coarse_states = spline(coarse)
    matches = []
    for index in range(0, len(times), 5):
        current = rows[index]
        def cost(phase, details=False):
            target_rows = spline(phase % period)
            squared = np.sum((current[:, None, :]-target_rows[None, :, :])**2, axis=2)
            ri, ci = linear_sum_assignment(squared)
            value = float(np.sum(squared[ri, ci]))
            return (math.sqrt(value), ci) if details else value
        values = []
        for q in coarse_states:
            sq = np.sum((current[:, None, :]-q[None, :, :])**2, axis=2)
            ri, ci = linear_sum_assignment(sq)
            values.append(sq[ri, ci].sum())
        guess = coarse[int(np.argmin(values))]
        fitted = minimize_scalar(cost, bounds=(guess-period/80, guess+period/80), method='bounded',
                                 options={'xatol': 1e-12})
        distance, permutation = cost(fitted.x, True)
        matches.append(dict(tau=float(times[index]), distance=distance, phase=float(fitted.x % period),
                            permutation_current_from_reference=permutation.tolist(),
                            contracting_pair=np.flatnonzero(np.isin(permutation, [2, 3])).tolist(),
                            expanding_pair=np.flatnonzero(np.isin(permutation, [6, 7])).tolist(),
                            halfcontrast_log10=np.log10(np.maximum(norms[index], np.nextafter(0., 1.))).tolist()))
    collapsed = np.all(rows[:, ::2, :] == rows[:, 1::2, :], axis=-1)
    recovered = contrast_norms(to_chart(data['full_states']))
    summary = json.loads((folder/'summary.json').read_text())
    recoveries, pending = [], []
    for pair in range(4):
        departures = next(e['times'] for e in summary['events']
                          if e['pair_index'] == pair and e['halfcontrast_threshold'] == 1e-4)
        low_indices = np.flatnonzero(norms[:, pair] < 1e-12)
        last_recovery = -1.
        for index in low_indices:
            if times[index] <= last_recovery:
                continue
            future = [t for t in departures if t > times[index]]
            if not future:
                pending.append(dict(pair=PAIRS[pair], first_low_checkpoint_tau=float(times[index]),
                                    minimum_halfcontrast=float(np.min(norms[index:, pair]))))
                break
            departure = future[0]
            mask = (times >= times[index]) & (times <= departure)
            recoveries.append(dict(pair=PAIRS[pair], first_low_checkpoint_tau=float(times[index]),
                                   minimum_halfcontrast=float(np.min(norms[mask, pair])),
                                   recovery_upward_crossing_1e_minus_4=departure))
            last_recovery = departure
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  minimum_halfcontrast_norms=norms.min(axis=0).tolist(),
                  minimum_log10_halfcontrast_norms=np.log10(np.maximum(norms.min(axis=0), np.nextafter(0., 1.))).tolist(),
                  final_log10_halfcontrast_norms=np.log10(np.maximum(norms[-1], np.nextafter(0., 1.))).tolist(),
                  nonzero_chart_at_every_checkpoint=bool(np.all(norms > 0)),
                  cartesian_reconstruction_equal_pair_checkpoints=np.sum(collapsed, axis=0).tolist(),
                  roundtrip_zero_contrast_checkpoints=np.sum(recovered == 0, axis=0).tolist(),
                  recovery_criterion='After a checkpoint below1e-12, first subsequent upward halfcontrast crossing of1e-4',
                  recovered_tiny_contrasts=recoveries, unresolved_small_contrasts_at_end=pending,
                  orbit_matches=matches, data_sha256=sha(folder/'trajectory.npz'), code_sha256=sha(__file__),
                  caveat='Near-zero saddle matching has interpolation and floating-point floors; chart contrasts are the precision diagnostic.')
    write(target, record)
    (folder/'analysis.source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps({k: record[k] for k in ['minimum_log10_halfcontrast_norms', 'final_log10_halfcontrast_norms',
                                            'nonzero_chart_at_every_checkpoint', 'roundtrip_zero_contrast_checkpoints']}))


def compare(args):
    first, second = args.output.resolve(), args.other.resolve()
    a, b = np.load(first/'trajectory.npz'), np.load(second/'trajectory.npz')
    ma = json.loads((first/'summary.json').read_text())
    mb = json.loads((second/'summary.json').read_text())
    exact_initial = np.array_equal(a['initial_chart_state'], b['initial_chart_state'])
    origin_a = a['origin_initial_chart_state'] if 'origin_initial_chart_state' in a else a['initial_chart_state']
    origin_b = b['origin_initial_chart_state'] if 'origin_initial_chart_state' in b else b['initial_chart_state']
    assert np.array_equal(origin_a, origin_b)
    assert np.array_equal(a['parameters'], b['parameters'])
    assert np.array_equal(a['time'], b['time'])
    events = []
    for e, f in zip(ma['events'], mb['events']):
        events.append(dict(pair=e['pair'], threshold=e['halfcontrast_threshold'],
                           first_tau=e['times'][0] if e['times'] else None,
                           second_tau=f['times'][0] if f['times'] else None,
                           first_event_difference=f['times'][0]-e['times'][0] if e['times'] and f['times'] else None,
                           number_of_events=[len(e['times']), len(f['times'])]))
    periods = np.array(a['time'])
    windows = []
    for end in [1000., 2500., 5000., 10000., 20000., 30000., 40000.]:
        mask = periods <= end
        if not np.any(mask) or end > periods[-1]:
            continue
        windows.append(dict(through_tau=end, maximum_chart_difference=float(np.max(abs(a['chart_states'][mask]-b['chart_states'][mask]))),
                            maximum_log10_norm_difference=float(np.max(abs(np.log10(a['halfcontrast_norms'][mask])-np.log10(b['halfcontrast_norms'][mask]))))))
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  folders=[first.relative_to(ROOT).as_posix(), second.relative_to(ROOT).as_posix()],
                  exactly_equal_initial_chart_states=bool(exact_initial), exactly_equal_origin_initial_chart_states=True,
                  exactly_equal_parameters=True, exactly_equal_checkpoint_times=True,
                  events=events, windows=windows, trajectory_hashes=[sha(first/'trajectory.npz'), sha(second/'trajectory.npz')],
                  source_sha256=sha(__file__), caveat='Long-time divergence is possible; event comparisons do not establish a heteroclinic theorem.')
    if (first/'diagnostics.json').exists() and (second/'diagnostics.json').exists():
        da = json.loads((first/'diagnostics.json').read_text())
        db = json.loads((second/'diagnostics.json').read_text())
        rec_a, rec_b = da['recovered_tiny_contrasts'], db['recovered_tiny_contrasts']
        recovery_comparison = []
        for ra in rec_a:
            candidates = [rb for rb in rec_b if rb['pair'] == ra['pair']]
            if not candidates:
                continue
            rb = min(candidates, key=lambda q: abs(q['recovery_upward_crossing_1e_minus_4']-ra['recovery_upward_crossing_1e_minus_4']))
            ta, tb = ra['recovery_upward_crossing_1e_minus_4'], rb['recovery_upward_crossing_1e_minus_4']
            recovery_comparison.append(dict(pair=ra['pair'], first_tau=ta, second_tau=tb, difference=tb-ta))
        record['recovery_event_comparison'] = recovery_comparison
    record['raw_event_caution'] = ('Repeated raw crossings, especially the0.01 threshold of already separated pairs, '
        'can be missed between accepted steps when a rapid oscillation dips narrowly below threshold. '
        'Their counts and first times are retained and are not used to classify saddle escapes. '
        'The recovery comparison requires a prior checkpoint below1e-12 followed by an upward1e-4 crossing.')
    write(first/'comparison.json', record)
    (first/'comparison.source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(record, indent=2))


def seed_compare(args):
    first, second = args.output.resolve(), args.other.resolve()
    a, b = np.load(first/'trajectory.npz'), np.load(second/'trajectory.npz')
    ma, mb = [json.loads((d/'summary.json').read_text()) for d in [first, second]]
    assert np.array_equal(a['time'][:len(b['time'])], b['time'])
    times = b['time']
    n1, n2 = a['halfcontrast_norms'][:len(times), 1], b['halfcontrast_norms'][:, 1]
    seed1, seed2 = a['halfcontrast_norms'][0, 1], b['halfcontrast_norms'][0, 1]
    offset = float(np.log10(seed1/seed2))
    period = ma['reference_period']
    strobes = abs(times/period-np.round(times/period)) < 1e-8
    fit_mask = strobes & (times >= 2500) & (times <= 5000)
    f1, f2 = np.polyfit(times[fit_mask], np.log(n1[fit_mask]), 1), np.polyfit(times[fit_mask], np.log(n2[fit_mask]), 1)
    growth = next(m['exponent'] for m in ma['initial_modes'] if m['pair'] == [6, 7])
    def first_departure(m):
        return next(e['times'][0] for e in m['events'] if e['pair'] == [2, 3] and e['halfcontrast_threshold'] == 1e-4)
    t1, t2 = first_departure(ma), first_departure(mb)
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  folders=[first.relative_to(ROOT).as_posix(), second.relative_to(ROOT).as_posix()],
                  seed_halfcontrast_norms=[float(seed1), float(seed2)], log10_seed_ratio=offset,
                  maximum_log10_ratio_error_through5000=float(np.max(abs(np.log10(n1[times <= 5000])-np.log10(n2[times <= 5000])-offset))),
                  stroboscopic_fit_window=[2500, 5000], fitted_growing_rates=[float(f1[0]), float(f2[0])],
                  reference_expanding_floquet_exponent=growth,
                  expected_logarithmic_delay=float(np.log(seed1/seed2)/growth),
                  first_upward_crossing_1e_minus_4=[t1, t2], measured_first_crossing_delay=t2-t1,
                  fitted_envelope_delay=float((f1[1]-f2[1])/((f1[0]+f2[0])/2)),
                  data_hashes=[sha(first/'trajectory.npz'), sha(second/'trajectory.npz')], code_sha256=sha(__file__),
                  caveat='Pointwise norm oscillates within the periodic orbit; a first threshold crossing is phase-dependent. The fitted envelope delay is a separate diagnostic.')
    write(first.parent/'seed-control.json', record)
    (first.parent/'seed-control.source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(record, indent=2))


def load_chain(folder):
    data = np.load(folder/'trajectory.npz')
    meta = json.loads((folder/'summary.json').read_text())
    times, states, norms = data['time'], data['chart_states'], data['halfcontrast_norms']
    matches = json.loads((folder/'diagnostics.json').read_text())['orbit_matches']
    if meta.get('resume'):
        pt, ps, pn, pm = load_chain(ROOT/meta['resume']['parent'])
        assert pt[-1] == times[0] and np.array_equal(ps[-1], states[0])
        times, states, norms = np.r_[pt, times[1:]], np.r_[ps, states[1:]], np.r_[pn, norms[1:]]
        matches = pm+matches[1:]
    return times, states, norms, matches


def plot(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder = args.output.resolve()
    times, states, norms, matches = load_chain(folder)
    fig, axes = plt.subplots(3 if args.seed_folder else 2, 1, figsize=(11, 9 if args.seed_folder else 7),
                             constrained_layout=True)
    colors = ['#4477AA', '#EE6677', '#228833', '#AA3377']
    for j, color in enumerate(colors):
        axes[0].plot(times, np.log10(norms[:, j]), color=color, lw=1.1,
                     label=f'Rows {2*j},{2*j+1}')
    axes[0].set(ylabel='log10 half-contrast norm', xlabel='Normalized time tau',
                title='Full 23-state chart: tiny row differences persist and recover')
    axes[0].legend(ncol=4, loc='lower left')
    axes[0].grid(alpha=.2)
    axes[1].semilogy([m['tau'] for m in matches], [m['distance'] for m in matches], color='#222222', lw=1)
    axes[1].axhline(1e-9, color='#888888', ls=':', label='Approximate interpolation floor')
    axes[1].set(ylabel='Distance to permuted saddle', xlabel='Normalized time tau',
                title='Finite-time visits to symmetry-related periodic saddles')
    axes[1].legend(loc='upper right')
    axes[1].grid(alpha=.2)
    if args.seed_folder:
        st, _, sn, _ = load_chain(args.seed_folder.resolve())
        if args.other:
            ct, _, cn, _ = load_chain(args.other.resolve())
        else:
            ct, cn = times, norms
        keep = ct <= st[-1]
        axes[2].plot(ct[keep], np.log10(cn[keep, 1]), label='Pair (2,3), initial half-contrast $10^{-12}$', color='#EE6677')
        axes[2].plot(st, np.log10(sn[:, 1]), label='Pair (2,3), initial half-contrast $10^{-25}$', color='#4477AA')
        axes[2].axhline(-4, color='#888888', ls=':', label='Recovery threshold $10^{-4}$')
        axes[2].set(xlabel='Normalized time tau', ylabel='log10 half-contrast norm', title='Smaller nonzero seeds delay the next departure')
        axes[2].legend(loc='lower right')
        axes[2].grid(alpha=.2)
    image = folder/'precision-diagnostic.png'
    fig.savefig(image, dpi=160)
    plt.close(fig)
    write(folder/'plot.json', dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                  image_sha256=sha(image), code_sha256=sha(__file__),
                                  data_sha256=sha(folder/'trajectory.npz'), plotted_time_end=float(times[-1])))
    (folder/'plot.source.py').write_bytes(Path(__file__).read_bytes())
    print(str(image.relative_to(ROOT)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['run', 'analyze', 'compare', 'seed-compare', 'plot'])
    parser.add_argument('--source', type=Path, default=ROOT/'results/motor-postcritical/six-cluster-cycle-rho-3p268/cycle.npz')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--other', type=Path)
    parser.add_argument('--resume-from', type=Path)
    parser.add_argument('--seed-folder', type=Path)
    parser.add_argument('--horizon', type=float, default=20000.)
    parser.add_argument('--chunk', type=float, default=500.)
    parser.add_argument('--max-step', type=float, default=.1)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--mean-atol', type=float, default=1e-12)
    parser.add_argument('--contrast-atol', type=float, default=1e-150)
    parser.add_argument('--unstable', type=float, default=1e-7)
    parser.add_argument('--seed', type=float, default=1e-12)
    parser.add_argument('--sign', type=int, choices=[-1, 1], default=1)
    args = parser.parse_args()
    if args.command == 'run':
        run(args)
    elif args.command == 'analyze':
        analyze(args)
    elif args.command == 'seed-compare':
        seed_compare(args)
    elif args.command == 'plot':
        plot(args)
    else:
        compare(args)

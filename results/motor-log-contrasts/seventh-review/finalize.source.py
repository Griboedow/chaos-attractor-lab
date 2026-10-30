"""Check and freeze the separate scalar continuation through tau280000."""
import datetime
from decimal import Decimal, localcontext
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parents[1]
sys.path.insert(0, str(ROOT/'research'))
import motor_log_contrasts as log
import motor_pair_coordinates as pair


def main():
    old_manifest = json.loads((BASE/'manifest-150k.json').read_text())
    for entry in old_manifest['files']:
        assert pair.sha(BASE/entry['path']) == entry['sha256']
    spec = importlib.util.spec_from_file_location('seventh_phase_review', HERE/'review.py')
    review = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(review)
    names = ['step-0p1-to-280k', 'step-0p05-to-280k']
    data, summaries, diagnostics, histories, validated = [], [], [], [], []
    for name in names:
        folder = BASE/name
        s = json.loads((folder/'summary.json').read_text())
        a = np.load(folder/'trajectory.npz')
        d = json.loads((folder/'diagnostics.json').read_text())
        parent = np.load(ROOT/s['parent']/'trajectory.npz')
        parent_summary = json.loads((ROOT/s['parent']/'summary.json').read_text())
        assert s['status'] == 'complete' and s['elapsed_tau'] == 280000.
        assert s['rhs_mode'] == 'scalar' and s['exact_saved_log_state_restored']
        assert a['initial_log_state'].tobytes() == parent['final_log_state'].tobytes()
        assert a['initial_log_state'].tobytes() == a['log_states'][0].tobytes()
        assert a['final_log_state'].tobytes() == a['log_states'][-1].tobytes()
        assert np.all(np.isfinite(a['log_states'])) and np.all(log.gauge_norm_squared(a['log_states']) > 0)
        for key, file in [('code_sha256', 'source.py'), ('pair_code_sha256', 'pair-source.py'),
                          ('scalar_code_sha256', 'scalar-source.py'), ('scalar_audit_sha256', 'rhs-audit.json'),
                          ('trajectory_sha256', 'trajectory.npz')]:
            assert s[key] == pair.sha(folder/file)
        assert s['parent_data_sha256'] == pair.sha(ROOT/s['parent']/'trajectory.npz')
        assert s['parent_summary_sha256'] == pair.sha(ROOT/s['parent']/'summary.json')
        assert d['data_sha256'] == pair.sha(folder/'trajectory.npz')
        history = review.history(folder)
        histories.append(history)
        times, chart, logs, _ = history
        minima = []
        for j in range(4):
            index = int(np.argmin(logs[:, j]))
            minima.append(dict(pair=pair.PAIRS[j], minimum_stored_log10_norm=float(logs[index, j]),
                               minimum_checkpoint_tau=float(times[index]), final_log10_norm=float(logs[-1, j])))
        rows = pair.chart_rows(chart)
        outside = []
        for i in range(8):
            for j in range(i+1, 8):
                if i//2 == j//2:
                    continue
                distance = np.linalg.norm(rows[:, i]-rows[:, j], axis=1)
                index = int(np.argmin(distance))
                outside.append(dict(rows=[i, j], minimum_sampled_distance=float(distance[index]), tau=float(times[index])))
        outside_min = min(outside, key=lambda item: item['minimum_sampled_distance'])
        fits = []
        t, normlogs = a['time'], a['log10_halfcontrast_norms']
        strobe = abs(t/s['reference_period']-np.round(t/s['reference_period'])) < 1e-8
        for left, right, growing, contracting in [(160000., 180000., 2, 0),
                                                   (220000., 240000., 2, 0),
                                                   (275000., 280000., 0, 3)]:
            mask = strobe & (t >= left) & (t <= right)
            centered = t[mask]-np.mean(t[mask])
            values = []
            for role, j in [('growing', growing), ('contracting', contracting)]:
                y = normlogs[mask, j]*log.LN10
                coefficient = np.polyfit(centered, y, 1)
                values.append(dict(role=role, pair=pair.PAIRS[j], rate_per_tau=float(coefficient[0]),
                                   maximum_log_amplitude_fit_residual=float(np.max(abs(y-np.polyval(coefficient, centered)))),
                                   minimum_log10_norm=float(np.min(normlogs[mask, j])),
                                   maximum_log10_norm=float(np.max(normlogs[mask, j]))))
            fits.append(dict(window=[left, right], fixed_phase_checkpoints=int(np.sum(mask)), fits=values,
                             measured_contraction_expansion_ratio=-values[1]['rate_per_tau']/values[0]['rate_per_tau']))
        bounds = []
        for item in s['underflow']:
            if item['maximum_omitted_component_log_bound'] is not None:
                with localcontext() as context:
                    context.prec = 60
                    bound = Decimal(str(item['maximum_omitted_component_log_bound'])).exp()
                    assert bound < Decimal('1e-322')
                bounds.append(dict(pair=item['pair'], maximum_recorded_log_bound=item['maximum_omitted_component_log_bound'],
                                   corresponding_decimal=str(bound)))
        validated.append(dict(run=name, status=s['status'], checkpoints=len(a['time']),
                              bitwise_exact27_state_resume=True, finite_logs_and_positive_gauge_norms=True,
                              wall_seconds=s['wall_seconds'], maximum_gauge_norm_squared_drift=s['maximum_gauge_norm_squared_drift'],
                              gauge_norm_squared_initial=log.gauge_norm_squared(a['initial_log_state']).tolist(),
                              gauge_norm_squared_final=log.gauge_norm_squared(a['final_log_state']).tolist(),
                              maximum_gauge_norm_squared_drift_since40k=max(parent_summary['maximum_gauge_norm_squared_drift'], float(np.max(abs(log.gauge_norm_squared(a['log_states'])-log.gauge_norm_squared(parent['initial_log_state']))))),
                              global_history_minima=minima, fixed_phase_rate_fits=fits,
                              sampled_nonfixed_pair_separation=outside_min,
                              omitted_mean_forcing_bounds=bounds,
                              recovered_tiny_contrasts=d['recovered_tiny_contrasts']))
        data.append(a)
        summaries.append(s)
        diagnostics.append(d)
    a, b = data
    assert np.array_equal(a['time'], b['time']) and np.array_equal(a['parameters'], b['parameters'])
    assert np.array_equal(a['origin_initial_chart_state'], b['origin_initial_chart_state'])
    windows = []
    for end in [175000., 200000., 225000., 250000., 275000., 280000.]:
        mask = a['time'] <= end
        za, zb = a['log_states'][mask], b['log_states'][mask]
        ua = np.stack([za[:, 11:15], za[:, 15:19], za[:, 19:23]], axis=-1)/np.sqrt(log.gauge_norm_squared(za))[..., None]
        ub = np.stack([zb[:, 11:15], zb[:, 15:19], zb[:, 19:23]], axis=-1)/np.sqrt(log.gauge_norm_squared(zb))[..., None]
        windows.append(dict(through_tau=end, maximum_mean_difference=float(np.max(abs(za[:, :11]-zb[:, :11]))),
                            maximum_log10_norm_difference=float(np.max(abs(a['log10_halfcontrast_norms'][mask]-b['log10_halfcontrast_norms'][mask]))),
                            maximum_normalized_direction_difference=float(np.max(abs(ua-ub))),
                            maximum_representable_chart_difference=float(np.max(abs(log.to_chart(za)-log.to_chart(zb))))))
    events = []
    for first in diagnostics[0]['recovered_tiny_contrasts']:
        candidates = [x for x in diagnostics[1]['recovered_tiny_contrasts'] if x['pair'] == first['pair']]
        second = min(candidates, key=lambda x: abs(x['recovery_upward_crossing_1e_minus_4']-first['recovery_upward_crossing_1e_minus_4']))
        ta, tb = first['recovery_upward_crossing_1e_minus_4'], second['recovery_upward_crossing_1e_minus_4']
        events.append(dict(pair=first['pair'], main_tau=ta, strict_tau=tb, difference=tb-ta))
    forecast = json.loads((BASE/'pre-seventh-prediction.json').read_text())
    comparison = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                      status='both280k_histories_verified', source_sha256=pair.sha(__file__),
                      frozen150k_manifest_sha256=pair.sha(BASE/'manifest-150k.json'), frozen150k_files_unchanged=True,
                      independent_pointwise_controls=True, phase_alignment_not_applied_to_control_comparison=True,
                      windows=windows, recovery_event_comparison=events, runs=validated,
                      unchanged_seventh_forecast=forecast,
                      forecast_residuals=[event['main_tau']-forecast['predicted_seventh_departure'] for event in events],
                      limitations=['Fixed-pair precision chart and finite-time evidence.',
                                   'Sampled row separation is not a continuous-time lower bound.',
                                   'Recorded omitted forcing bounds do not bound accumulated trajectory error.',
                                   'No heteroclinic-connection, attracting-network, or chaos theorem.'])
    pair.write(HERE/'comparison.json', comparison)
    (HERE/'finalize.source.py').write_bytes(Path(__file__).read_bytes())
    paths = []
    for name in names+['seventh-review']:
        paths.extend(p for p in (BASE/name).rglob('*') if p.is_file() and p.name != 'manifest.json' and '__pycache__' not in p.parts)
    paths.append(BASE/'seventh-progress.json')
    entries = [dict(path=p.relative_to(BASE).as_posix(), bytes=p.stat().st_size, sha256=pair.sha(p)) for p in sorted(paths)]
    pair.write(HERE/'manifest.json', dict(created_utc=comparison['created_utc'], status=comparison['status'],
                                         scope='Separate280k continuation package; does not replace150k evidence.',
                                         frozen150k_manifest_sha256=comparison['frozen150k_manifest_sha256'],
                                         forecast_sha256=pair.sha(BASE/'pre-seventh-prediction.json'),
                                         total_bytes=sum(p['bytes'] for p in entries), file_count=len(entries), files=entries))
    print(json.dumps(dict(status=comparison['status'], windows=windows, events=events,
                          minima=[r['global_history_minima'] for r in validated], file_count=len(entries)), indent=2))


if __name__ == '__main__':
    main()

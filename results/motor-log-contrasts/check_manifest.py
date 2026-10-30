"""Freeze only the two completed150k histories and their supporting audits."""
import datetime
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import scipy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    scopes = ['step-0p1', 'step-0p05', 'scalar-benchmark', 'scalar-benchmark-recheck', 'orbit-matching-review']
    runs, data, bounds = [], [], []
    for name in scopes[:2]:
        folder = HERE/name
        record = json.loads((folder/'summary.json').read_text())
        assert record['status'] == 'complete'
        assert record['elapsed_tau'] == 150000.
        assert record['code_sha256'] == sha(folder/'source.py')
        assert record['pair_code_sha256'] == sha(folder/'pair-source.py')
        assert record['trajectory_sha256'] == sha(folder/'trajectory.npz')
        parent = ROOT/record['parent']
        assert record['parent_data_sha256'] == sha(parent/'trajectory.npz')
        assert record['parent_summary_sha256'] == sha(parent/'summary.json')
        a, p = np.load(folder/'trajectory.npz'), np.load(parent/'trajectory.npz')
        assert np.array_equal(a['initial_chart_state'], p['final_chart_state'])
        assert np.array_equal(a['initial_log_state'], a['log_states'][0])
        assert np.array_equal(a['final_log_state'], a['log_states'][-1])
        assert np.all(np.diff(a['time']) > 0)
        assert np.all(np.isfinite(a['log_states']))
        for item in record['underflow']:
            if item['maximum_omitted_component_log_bound'] is not None:
                with localcontext() as context:
                    context.prec = 60
                    converted = Decimal(str(item['maximum_omitted_component_log_bound'])).exp()
                    assert converted < Decimal('1e-322')
                bounds.append(dict(run=name, pair=item['pair'],
                                   maximum_recorded_log_bound=item['maximum_omitted_component_log_bound'],
                                   corresponding_decimal=str(converted),
                                   scope='Instantaneous algebraic forcing bound at recorded floating-point stages; not an interval trajectory-error certificate.'))
        diagnostic = json.loads((folder/'diagnostics.json').read_text())
        assert diagnostic['data_sha256'] == sha(folder/'trajectory.npz')
        assert diagnostic['finite_log_states'] and diagnostic['positive_gauge_norms']
        runs.append(dict(name=name, status=record['status'], final_tau=150000., checkpoints=len(a['time']),
                         maximum_gauge_norm_squared_drift=record['maximum_gauge_norm_squared_drift'],
                         finite_log_states=True, source_hashes_verified=True, initial_saved_chart_verified=True))
        data.append(a)
    assert np.array_equal(data[0]['time'], data[1]['time'])
    assert np.array_equal(data[0]['parameters'], data[1]['parameters'])
    assert np.array_equal(data[0]['origin_initial_chart_state'], data[1]['origin_initial_chart_state'])
    for name in scopes[2:4]:
        folder = HERE/name
        audit = json.loads((folder/'audit.json').read_text())
        assert audit['status'] == 'passed'
        assert audit['code_sha256'] == sha(folder/'source.py')
        assert audit['log_code_sha256'] == sha(folder/'log-source.py')
        assert audit['pair_code_sha256'] == sha(folder/'pair-source.py')
    review = json.loads((HERE/'orbit-matching-review/audit.json').read_text())
    assert review['status'] == 'completed'
    assert review['source_sha256'] == sha(HERE/'orbit-matching-review/source.py')
    assert review['numerical_trajectory_unchanged']
    assert review['trajectory_sha256'] == sha(HERE/'step-0p1/trajectory.npz')
    original = json.loads((HERE/'pre-event-prediction.json').read_text())
    extension = json.loads((HERE/'pre-seventh-prediction.json').read_text())
    assert extension['source_sha256'] == sha(HERE/'pre-event-prediction.json')
    next_gap = original['kappa']*(original['predicted_sixth_departure']-original['predicted_fifth_departure'])+original['offset_from_last_two_gaps']
    assert extension['predicted_seventh_departure'] == original['predicted_sixth_departure']+next_gap
    names = ['pre-event-prediction.json', 'pre-seventh-prediction.json', 'frozen-40k-node-visits.json',
             'frozen-seed-1e-25-node-visits.json', 'check_manifest.py', 'review_orbit_matching.py']
    files = [HERE/name for name in names]
    for name in scopes:
        files.extend(p for p in (HERE/name).rglob('*') if p.is_file())
    entries = [dict(path=p.relative_to(HERE).as_posix(), bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(files)]
    result = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  status='completed150k_histories_verified', scope_directories=scopes,
                  excludes_later_continuation_folders=True,
                  environment=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__),
                  runs=runs, identical_checkpoint_times_parameters_and_origin=True,
                  observed_omitted_mean_forcing_bounds=bounds,
                  total_bytes=sum(e['bytes'] for e in entries), file_count=len(entries), files=entries)
    (HERE/'manifest-150k.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: result[k] for k in ['status', 'total_bytes', 'file_count', 'runs']}, indent=2))


if __name__ == '__main__':
    main()

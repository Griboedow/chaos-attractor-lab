"""Freeze the completed CPU evidence inventory, leaving the shared report editable."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    folder = Path(__file__).resolve().parent
    out = folder / 'manifest.json'
    if out.exists():
        raise FileExistsError('Do not overwrite a frozen inventory')
    run = folder / 'run-01'
    record = json.loads((run / 'analysis.json').read_text())
    derived = json.loads((folder / 'analysis/analysis.json').read_text())
    strict = json.loads((folder / 'strict-control.json').read_text())
    old = json.loads((run / 'input.json').read_text())
    checks = dict(completed=record['status'] == 'complete', observed_exactly_20000=record['observed_tau'] == 20000.,
                  added_exactly_16000=record['added_tau'] == 16000., sixteen_complete_blocks=len(record['new_blocks']) == 16 and all(b['complete'] and b['duration'] == 1000. for b in record['new_blocks']),
                  endpoint_controls=record['initial_endpoint_control']['passed'] and record['final_endpoint_control']['passed'],
                  independent_audit=derived['status'] == 'passed' and all(derived['controls'].values()),
                  strict_control=strict['status'] == 'passed' and all(strict['controls'].values()),
                  strict_source=sha(folder / 'strict-control.source.py') == strict['source_sha256'],
                  analysis_source=sha(folder / 'analyze.source.py') == derived['source_sha256'],
                  analysis_record=sha(run / 'analysis.json') == derived['input_report_sha256'],
                  unchanged_input_files=record['source_inputs_unchanged'],
                  inherited_windows=record['windows'][:len(old['windows'])] == old['windows'])
    for name, digest in strict['input_hashes'].items():
        checks['strict_input_' + name] = sha(run / name) == digest
    with np.load(run / 'input.npz') as initial, np.load(run / 'continuation.npz') as final:
        previous_time = float(initial['elapsed_original'])
        restored_logs = initial['log_sums'].copy()
        restored_trace = float(initial['integrated_trace'])
        speed = 0.
        for block in record['new_blocks']:
            name = f"block-{block['index']:02d}"
            block_json = json.loads((run / (name + '.json')).read_text())
            checks[name + '_metadata'] = block_json == block
            checks[name + '_npz_hash'] = sha(run / (name + '.npz')) == block['npz_sha256']
            with np.load(run / (name + '.npz')) as part:
                times, states = part['time'], part['states']
                restored_logs += part['block_log_sums']
                restored_trace += float(part['block_integrated_trace'])
                speed += float(part['block_integrated_speed_squared'])
                checks[name + '_grid'] = bool(len(times) == 5000 and abs(times[0] - previous_time - .2) < 1e-10 and np.max(abs(np.diff(times) - .2)) < 1e-10)
                checks[name + '_state_endpoint'] = bool(np.array_equal(states[-1], part['final_state']))
                previous_time = float(times[-1])
        # Re-grouping changes floating-point summation; retain measured discrepancies.
        discrepancies = dict(log_sums=float(np.max(abs(restored_logs - final['log_sums']))),
                             integrated_trace=abs(restored_trace - float(final['integrated_trace'])),
                             new_speed_integral=abs(speed - float(final['integrated_speed_squared_new'])))
        checks['grouped_accumulators_agree'] = max(discrepancies.values()) < 1e-7
        checks['final_time_matches'] = previous_time == float(final['elapsed_original'])
        counts = {name: list(final[name].shape) for name in ('time', 'states', 'section_time', 'section_states', 'checkpoint')}
        checks['sample_count'] = counts['time'] == [101000] and counts['states'] == [101000, 23]
        checks['checkpoint_shape'] = counts['checkpoint'] == [58, 579]
    assert all(checks.values()), {name: value for name, value in checks.items() if not value}
    files = {path.relative_to(folder).as_posix(): dict(bytes=path.stat().st_size, sha256=sha(path))
             for path in sorted(folder.rglob('*')) if path.is_file() and '__pycache__' not in path.parts}
    result = dict(status='passed', created_utc=datetime.now(timezone.utc).isoformat(), checks=checks,
                  grouped_accumulator_discrepancies=discrepancies, history_shapes=counts,
                  files=files, file_count=len(files), total_bytes=sum(v['bytes'] for v in files.values()),
                  scope='Immutable CPU evidence only. The shared docs/motor-activity-continuation.md may receive independently owned GPU and cluster sections and is deliberately outside this inventory.')
    out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({key: result[key] for key in ('status', 'file_count', 'total_bytes', 'history_shapes')}))


if __name__ == '__main__':
    main()

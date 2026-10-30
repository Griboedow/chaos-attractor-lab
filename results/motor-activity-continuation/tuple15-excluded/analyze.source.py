"""Audit both retained cases and compare them with the fitted periodic saddle."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment, minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FOLDER = Path(__file__).resolve().parent
ROOT = FOLDER.parents[2]
OUT = FOLDER / 'analysis'
if OUT.exists():
    raise FileExistsError('Preserve previous analysis')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
spec = importlib.util.spec_from_file_location('frozen_rows', FOLDER / 'lane-47/motor_rows.py')
motor = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = motor
spec.loader.exec_module(motor)
records, samples, checks = {}, {}, {}
for lane in (45, 47):
    run = FOLDER / f'lane-{lane}'
    record = json.loads((run / 'analysis.json').read_text())
    assert 'periodic_attempts' in record
    records[lane] = record
    checks[f'{lane}_complete'] = record['status'] == 'complete' and record['elapsed_added_tau'] == 500.
    checks[f'{lane}_clock'] = record['refined_segment_elapsed'] == 1500. and record['scout_elapsed'] == 600. and record['cumulative_source_elapsed'] == 2100. and record['cumulative_elapsed_tau'] == 2600.
    checks[f'{lane}_sample_hash'] = sha(run / 'sample.npz') == record['sample_sha256']
    checks[f'{lane}_seed_hash'] = sha(run / 'seed.npz') == record['seed_sha256']
    checks[f'{lane}_source_hash'] = sha(ROOT / record['source_path']) == record['input_source_sha256']
    checks[f'{lane}_scout_hash'] = sha(ROOT / record['scout_path']) == record['scout_source_sha256']
    for name, digest in record['source_hashes'].items():
        checks[f'{lane}_source_{name}'] = sha(run / name) == digest
    with np.load(run / 'sample.npz') as data:
        samples[lane] = {name: data[name].copy() for name in data.files}
    data = samples[lane]
    checks[f'{lane}_grid'] = bool(np.array_equal(data['time'], np.arange(1001) / 2))
    checks[f'{lane}_no_state_reset'] = bool(np.array_equal(data['states'][0], record['initial_state']))
    checks[f'{lane}_basis'] = bool(np.max(abs(data['final_basis'].T @ data['final_basis'] - np.eye(23))) < 1e-12)
    checks[f'{lane}_trace'] = record['trace_error_per_tau'] < 1e-9 and record['max_local_trace_error'] < 1e-8
    checks[f'{lane}_flow'] = record['max_local_flow_error'] < 1e-8 and abs(record['flow_exponent'] - record['telescopic_flow_exponent']) < 1e-9
    checks[f'{lane}_endpoint'] = max(record['final_endpoint_control'].values()) < 1e-10
    checks[f'{lane}_all_lags'] = [r['lag'] for r in record['return_lags']] == list(range(1, 17))
periodic = json.loads((FOLDER / 'lane-47/periodic/analysis.json').read_text())
checks['periodic_complete'] = periodic['status'] == 'complete'
checks['periodic_closures'] = max(periodic[k] for k in ('closure_error', 'Radau_closure_error', 'full_closure_error')) < 1e-9
checks['periodic_flow'] = periodic['full_flow_identity_error'] < 1e-8
checks['periodic_independent_variation'] = periodic['DOP853_Radau_variational_control']['error'] < 1e-8
with np.load(FOLDER / 'lane-47/periodic/cycle.npz') as data:
    p, initial, period = data['parameters'], data['full_states'][0], float(data['time'][-1])
    multipliers = np.linalg.eigvals(data['full_monodromy'])
    quotient = np.linalg.eigvals(data['quotient_monodromy'])
checks['three_unstable_full_multipliers'] = int(np.sum(abs(multipliers) > 1 + 1e-6)) == 3
checks['one_unit_full_multiplier'] = int(np.sum(abs(abs(multipliers) - 1) < 1e-6)) == 1
quotient_nonflow = np.delete(quotient, np.argmin(abs(quotient - 1)))
checks['stable_quotient_nonflow'] = bool(np.max(abs(quotient_nonflow)) < 1)
solution = solve_ivp(lambda t, y: motor.rhs(t, y, p, 8), (0, period), initial, method='DOP853',
                     rtol=2e-13, atol=2e-15, max_step=.01, dense_output=True)
checks['reference_flow_success'] = solution.success
checks['reference_closure'] = bool(np.max(abs(solution.y[:, -1] - initial)) < 1e-9)


def physical(state):
    return np.column_stack((state[:8], state[8:16], np.r_[state[16:], -state[16:].sum()]))


def distance(current, phase):
    reference = physical(solution.sol(phase % period))
    cost = np.sum((current[:, None, :] - reference[None, :, :]) ** 2, axis=2)
    a, b = linear_sum_assignment(cost)
    permutation = np.empty(8, dtype=int)
    permutation[a] = b
    return float(np.sqrt(cost[a, b].sum())), permutation.tolist()


grid = np.arange(2000) * period / 2000
reference = np.array([physical(solution.sol(t)) for t in grid])
step = period / 2000
fits = {}
for lane in (45, 47):
    fits[lane] = []
    for index in range(0, 1001, 50):
        state = samples[lane]['states'][index]
        current = physical(state)
        costs = np.sum((current[None, :, None, :] - reference[:, None, :, :]) ** 2, axis=3)
        rough = []
        for cost in costs:
            a, b = linear_sum_assignment(cost)
            rough.append(np.sqrt(cost[a, b].sum()))
        rough = np.array(rough)
        minima = np.flatnonzero((rough <= np.roll(rough, 1)) & (rough <= np.roll(rough, -1)))
        candidates = sorted(minima, key=lambda j: rough[j])[:8]
        best = None
        for j in candidates:
            result = minimize_scalar(lambda t: distance(current, t)[0], bounds=(grid[j] - step, grid[j] + step),
                                     method='bounded', options=dict(xatol=1e-13, maxiter=120))
            value, permutation = distance(current, result.x)
            entry = dict(added_tau=float(samples[lane]['time'][index]), distance=value,
                         phase=float(result.x % period), current_from_reference_permutation=permutation)
            if best is None or entry['distance'] < best['distance']:
                best = entry
        fits[lane].append(best)
assert all(checks.values()), {key: value for key, value in checks.items() if not value}
OUT.mkdir()
fig, ax = plt.subplots(figsize=(8, 4.5))
for lane, color in ((45, '#A85636'), (47, '#477D99')):
    ax.semilogy([v['added_tau'] for v in fits[lane]], [v['distance'] for v in fits[lane]], 'o-',
                color=color, label=f'Source lane {lane}')
ax.set(xlabel='Added original tau; no warmup', ylabel='Distance to saddle orbit, allowing row permutation',
       title='Same parameter: near-saddle escape and periodic capture')
ax.grid(alpha=.2)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / 'orbit-distance.png', dpi=170)
plt.close(fig)
result = dict(status='passed', created_utc=datetime.now(timezone.utc).isoformat(), checks=checks,
              source_sha256=sha(Path(__file__)), input_analysis_hashes={str(lane): sha(FOLDER / f'lane-{lane}/analysis.json') for lane in (45, 47)},
              periodic_analysis_sha256=sha(FOLDER / 'lane-47/periodic/analysis.json'),
              periodic_cycle_sha256=sha(FOLDER / 'lane-47/periodic/cycle.npz'),
              period=period, unstable_exponents=sorted(np.log(abs(multipliers[abs(multipliers) > 1 + 1e-6])) / period, reverse=True),
              quotient_max_nonflow=float(np.max(abs(quotient_nonflow))), phase_permutation_fits=fits,
              phase_matching=dict(grid_points=2000, refined_local_minima=8, reference_max_step=.01,
                                  reference_rtol=2e-13, reference_atol=2e-15,
                                  convention='current row j approximates reference row permutation[j] at the stored phase. Distances use all 24 physical row coordinates.'),
              eigenvalue_limitation='Strongly contracting multipliers over period62.77 approach floating-point noise. Their individual negative exponents are not accurately resolved by this monodromy calculation. The three expanding multipliers are well separated from the unit circle.',
              figure_sha256=sha(OUT / 'orbit-distance.png'),
              limitation='Finite numerical classification of two archived states; no asymptotic chaos, attraction or rigorous Floquet theorem.')
(OUT / 'analysis.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
inventory = {path.relative_to(FOLDER).as_posix(): dict(bytes=path.stat().st_size, sha256=sha(path))
             for path in sorted(FOLDER.rglob('*')) if path.is_file() and '__pycache__' not in path.parts}
(FOLDER / 'manifest.json').write_text(json.dumps(dict(status='passed', checks=checks, files=inventory,
    file_count=len(inventory)), indent=2) + '\n', encoding='utf-8', newline='\n')
print(json.dumps(dict(status='passed', period=period, unstable_exponents=result['unstable_exponents'],
    initial_fit_distances={lane: fit[0]['distance'] for lane, fit in fits.items()},
    final_fit_distances={lane: fit[-1]['distance'] for lane, fit in fits.items()})))

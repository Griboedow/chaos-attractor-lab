"""Audit and freeze the complete periodic-branch continuation and refinement."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
out = root / 'audit'
if out.exists():
    raise FileExistsError('Preserve previous branch audit')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
branch = json.loads((root / 'analysis.json').read_text())
refinement = json.loads((root / 'crossing-refinement/analysis.json').read_text())
checks = dict(branch_complete=branch['status'] == 'complete', target_reached=branch['final_successful_rho'] == branch['target_rho'],
              refinement_complete=refinement['status'] == 'complete',
              refinement_source_hash=sha(root / 'refine.source.py') == refinement['source_sha256'],
              refinement_parent_hash=sha(root / 'analysis.json') == refinement['branch_analysis_sha256'],
              input_cycle_hash=sha(root / 'input-cycle.npz') == branch['input_cycle_sha256'],
              input_metadata_hash=sha(root / 'input-analysis.json') == branch['input_analysis_sha256'])
for name, digest in branch['source_hashes'].items():
    checks['source_' + name] = sha(root / name) == digest
points = []
previous = None
for point in branch['attempts']:
    run = root / f"attempt-{point['index']:02d}"
    checks[f"attempt_{point['index']}_record"] = json.loads((run / 'analysis.json').read_text()) == point
    if previous is not None:
        checks[f"attempt_{point['index']}_previous_state"] = point['input_state'] == previous['initial_state'] and point['input_period'] == previous['period']
    if point['status'] != 'complete':
        continue
    points.append(point)
    checks[f"point_{point['index']}_hash"] = sha(run / 'cycle.npz') == point['cycle_sha256']
    with np.load(run / 'cycle.npz') as data:
        checks[f"point_{point['index']}_parameters"] = np.array_equal(data['parameters'], [point['rho'], *branch['fixed_r_b']])
        checks[f"point_{point['index']}_weights"] = np.array_equal(data['weights'], branch['weights'])
        checks[f"point_{point['index']}_state"] = np.array_equal(data['initial_state'], point['initial_state'])
        q = np.linalg.eigvals(data['quotient_monodromy'])
        q = np.delete(q, np.argmin(abs(q - 1)))
        c = [np.linalg.eigvals(m) for m in data['contrast_monodromies']]
        checks[f"point_{point['index']}_quotient_spectrum"] = abs(np.max(abs(q)) - point['spectrum']['quotient_max_nonflow']) < 1e-9
        checks[f"point_{point['index']}_contrast_spectra"] = np.max(abs(np.array([np.max(abs(v)) for v in c]) - point['spectrum']['contrast_maxima'])) < 1e-9
    checks[f"point_{point['index']}_controls"] = max(point['spectrum']['errors'].values()) < 1e-7 and point['strict_shooting_matrix_difference'] < 1e-7
    previous = point
checks['descending_successful_parameters'] = bool(np.all(np.diff([p['rho'] for p in points]) < 0))
checks['eighteen_successful_points'] = len(points) == 18
checks['no_transverse_instability_observed'] = all(p['spectrum']['contrast_unstable_count'] == 0 for p in points)
checks['two_quotient_unstable_at_target'] = points[-1]['spectrum']['quotient_unstable_count'] == 2
for call in refinement['calls']:
    path = root / 'crossing-refinement' / call['npz_path']
    checks[f"refinement_{call['index']}_hash"] = sha(path) == call['npz_sha256']
    checks[f"refinement_{call['index']}_controls"] = call['closure'] < 1e-8 and call['strict_matrix_difference'] < 1e-7
    with np.load(path) as data:
        values = np.linalg.eigvals(data['quotient_monodromy'])
        values = np.delete(values, np.argmin(abs(values - 1)))
        checks[f"refinement_{call['index']}_residual"] = abs(np.max(abs(values)) - 1 - call['residual']) < 1e-9
bracket = refinement['numerical_sign_bracket']
checks['numerical_sign_change'] = bracket['lower_value'] > 0 and bracket['upper_value'] < 0
checks['critical_complex_pair'] = abs(refinement['multiplier_at_threshold']['imag']) > .1 and abs(refinement['multiplier_at_threshold']['modulus'] - 1) < 1e-8
checks['sixth_period_symmetry'] = all(p['maximum_component_residual'] < 1e-9 for p in refinement['sixth_period_symmetry_controls'])
checks = {k: bool(v) for k, v in checks.items()}
assert all(checks.values()), {k: v for k, v in checks.items() if not v}
rho = np.array([p['rho'] for p in points])
qmax = np.array([p['spectrum']['quotient_max_nonflow'] for p in points])
cmax = np.array([max(p['spectrum']['contrast_maxima']) for p in points])
critical = refinement['numerical_threshold']
fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5))
axes[0].semilogy(rho, qmax, 'o-', label='Largest quotient nonflow modulus', color='#A85636')
axes[0].semilogy(rho, cmax, 'o-', label='Largest of all 9 pair-contrast moduli', color='#477D99')
axes[0].axhline(1, color='black', linewidth=.8)
axes[0].axvline(critical, color='#777777', linestyle='--', linewidth=.9)
axes[0].set(xlabel='rho', ylabel='Floquet multiplier modulus', title='Complete periodic branch: fixed cluster weights')
axes[0].legend(fontsize=8, loc='best')
axes[0].grid(alpha=.2)
theta = np.linspace(0, 2*np.pi, 401)
axes[1].plot(np.cos(theta), np.sin(theta), color='#777777', linestyle='--', label='Unit circle')
selected = [p['spectrum']['quotient_nonflow_multipliers'][0] for p in points if p['spectrum']['quotient_max_nonflow'] <= 1.5]
real = np.array([v['real'] for v in selected])
imag = np.abs([v['imag'] for v in selected])
axes[1].plot(real, imag, 'o-', color='#A85636', label='Leading quotient pair')
axes[1].plot(real, -imag, 'o-', color='#A85636')
mu = refinement['multiplier_at_threshold']
axes[1].scatter([mu['real']]*2, [mu['imag'], -mu['imag']], marker='*', s=110, color='#1B2028', label='Numerical crossing')
axes[1].set(xlabel='Real part', ylabel='Imaginary part', title='Complex pair near crossing (modulus <= 1.5)', aspect='equal')
axes[1].legend(fontsize=8, loc='lower right')
axes[1].grid(alpha=.2)
fig.suptitle('Five-cluster cycle: quotient instability while pair contrasts contract', fontsize=12)
fig.tight_layout()
out.mkdir()
fig.savefig(out / 'branch.png', dpi=170)
plt.close(fig)
summary = dict(status='passed', created_utc=datetime.now(timezone.utc).isoformat(), checks=checks,
               source_sha256=sha(Path(__file__)), branch_analysis_sha256=sha(root / 'analysis.json'),
               refinement_analysis_sha256=sha(root / 'crossing-refinement/analysis.json'), successful_points=len(points),
               failed_attempts=len(branch['attempts'])-len(points), numerical_threshold=critical,
               critical_period=refinement['period_at_threshold'], critical_multiplier=refinement['multiplier_at_threshold'],
               target=dict(rho=points[-1]['rho'], period=points[-1]['period'], spectrum=points[-1]['spectrum']),
               maximum_closure_error=max(p['spectrum']['errors']['closure'] for p in points),
               maximum_strict_matrix_difference=max(p['strict_shooting_matrix_difference'] for p in points),
               figure_sha256=sha(out / 'branch.png'),
               limitation='Numerical periodic-branch and complex multiplier crossing evidence, not a bifurcation theorem or a demonstrated connection to the chaotic candidate.')
(out / 'analysis.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8', newline='\n')
inventory = {p.relative_to(root).as_posix(): dict(bytes=p.stat().st_size, sha256=sha(p))
             for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
(root / 'manifest.json').write_text(json.dumps(dict(status='passed', checks=checks, files=inventory,
    file_count=len(inventory)), indent=2) + '\n', encoding='utf-8', newline='\n')
print(json.dumps(dict(status='passed', checks=len(checks), files=len(inventory), threshold=critical,
    maximum_closure=summary['maximum_closure_error'], maximum_matrix_difference=summary['maximum_strict_matrix_difference'])))

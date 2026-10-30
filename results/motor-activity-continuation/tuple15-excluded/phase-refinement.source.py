"""Refine the scalar phase after bounded search; preserve its original output."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
out = root / 'phase-refinement'
if out.exists():
    raise FileExistsError('Preserve previous refinement')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
original = json.loads((root / 'analysis/analysis.json').read_text())
spec = importlib.util.spec_from_file_location('frozen_rows_phase', root / 'lane-47/motor_rows.py')
motor = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = motor
spec.loader.exec_module(motor)
with np.load(root / 'lane-47/periodic/cycle.npz') as data:
    p, initial, period = data['parameters'], data['full_states'][0], float(data['time'][-1])
solution = solve_ivp(lambda t, y: motor.rhs(t, y, p, 8), (0, period), initial, method='DOP853',
                     rtol=2e-13, atol=2e-15, max_step=.01, dense_output=True)
assert solution.success
def physical(z):
    return np.column_stack((z[:8], z[8:16], np.r_[z[16:], -z[16:].sum()]))
fits = {}
for lane in (45, 47):
    with np.load(root / f'lane-{lane}/sample.npz') as data:
        states = data['states'].copy()
    fits[lane] = []
    for old in original['phase_permutation_fits'][str(lane)]:
        current = physical(states[int(2 * old['added_tau'])])
        phase, permutation = old['phase'], old['current_from_reference_permutation']
        for _ in range(12):
            state = solution.sol(phase % period)
            delta = current - physical(state)[permutation]
            velocity = physical(motor.rhs(0, state, p, 8))[permutation]
            correction = float(np.sum(delta * velocity) / np.sum(velocity * velocity))
            candidate = (phase + correction) % period
            if abs(correction) > .001 or np.linalg.norm(current - physical(solution.sol(candidate))[permutation]) > np.linalg.norm(delta):
                break
            phase = candidate
            if abs(correction) < 1e-14:
                break
        distance = float(np.linalg.norm(current - physical(solution.sol(phase))[permutation]))
        assert distance <= old['distance'] + 1e-12
        fits[lane].append(dict(added_tau=old['added_tau'], distance=distance, phase=float(phase),
                              current_from_reference_permutation=permutation))
out.mkdir()
fig, ax = plt.subplots(figsize=(8, 4.5))
for lane, color in ((45, '#A85636'), (47, '#477D99')):
    ax.semilogy([v['added_tau'] for v in fits[lane]], [v['distance'] for v in fits[lane]], 'o-', color=color, label=f'Source lane {lane}')
ax.set(xlabel='Added original tau; no warmup', ylabel='Distance to saddle orbit, allowing row permutation',
       title='Same parameter: near-saddle escape and periodic capture')
ax.grid(alpha=.2)
ax.legend()
fig.tight_layout()
fig.savefig(out / 'orbit-distance.png', dpi=170)
plt.close(fig)
early = [v for v in fits[45] if v['added_tau'] <= 150]
slope, intercept = np.polyfit([v['added_tau'] for v in early], np.log([v['distance'] for v in early]), 1)
record = dict(status='passed', created_utc=datetime.now(timezone.utc).isoformat(), source_sha256=sha(Path(__file__)),
              input_analysis_sha256=sha(root / 'analysis/analysis.json'), phase_permutation_fits=fits,
              exploratory_early_escape_log_distance_slope=float(slope),
              slope_scope='Post-hoc descriptive straight-line fit to log distance at added tau0,25,...150; periodic phase modulation is not modeled. Not an independent theorem or prospective prediction.',
              method='At most12 safeguarded scalar Gauss-Newton steps from each frozen bounded-search phase, holding its optimal row assignment. Reject steps above0.001 tau or increasing the residual.',
              figure_sha256=sha(out / 'orbit-distance.png'))
(out / 'analysis.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8', newline='\n')
inventory = {p.relative_to(root).as_posix(): dict(bytes=p.stat().st_size, sha256=sha(p))
             for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
(root / 'manifest-final.json').write_text(json.dumps(dict(status='passed', files=inventory,
    file_count=len(inventory), checks=original['checks']), indent=2) + '\n', encoding='utf-8', newline='\n')
print(json.dumps(dict(initial_distances={lane: values[0]['distance'] for lane, values in fits.items()},
                     final_distances={lane: values[-1]['distance'] for lane, values in fits.items()}, early_slope=float(slope))))

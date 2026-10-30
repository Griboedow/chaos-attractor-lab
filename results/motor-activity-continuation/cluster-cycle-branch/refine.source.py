"""Refine the numerical complex-pair threshold on the saved periodic branch."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import time
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.optimize import brentq

root = Path(__file__).resolve().parent
out = root / 'crossing-refinement'
if out.exists():
    raise FileExistsError('Preserve previous refinement')
sys.path.insert(0, str(root))
import motor_postcritical_cycles as cycles
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
branch = json.loads((root / 'analysis.json').read_text())
assert branch['status'] == 'complete'
points = [p for p in branch['attempts'] if p['status'] == 'complete']
brackets = [(a, b) for a, b in zip(points[:-1], points[1:])
            if (a['spectrum']['quotient_max_nonflow'] - 1) * (b['spectrum']['quotient_max_nonflow'] - 1) < 0]
assert len(brackets) == 1
above_rho, below_rho = brackets[0]
weights = np.array(branch['weights'])
groups = branch['groups']
seeds = [(p['rho'], np.array(p['initial_state']), p['period']) for p in points]
calls = []
out.mkdir()
started = time.perf_counter()
record = dict(status='running', started_utc=datetime.now(timezone.utc).isoformat(),
              branch_analysis_sha256=sha(root / 'analysis.json'), source_sha256=sha(Path(__file__)),
              bracket=[below_rho['rho'], above_rho['rho']], calls=calls,
              limitation='Floating-point root refinement of a complex quotient multiplier modulus. No interval enclosure, nondegenerate bifurcation theorem or invariant torus is proved.')


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8', newline='\n')


def evaluate(rho):
    closest = min(seeds, key=lambda p: abs(p[0] - rho))
    p = np.array([rho, *branch['fixed_r_b']])
    state, period, shooting, history = cycles.shoot(closest[1], closest[2], p, weights)
    solution = cycles.integrate(state, period, p, weights, augmented=True, rtol=2e-13, max_step=.01)
    d = len(state)
    matrix = solution.y[d:, -1].reshape(d, d)
    values = np.linalg.eigvals(matrix)
    flow_index = int(np.argmin(abs(values - 1)))
    nonflow = np.delete(values, flow_index)
    dominant = nonflow[np.argmax(abs(nonflow))]
    if dominant.imag < 0:
        dominant = dominant.conjugate()
    assert abs(dominant.imag) > .1
    residual = float(abs(dominant) - 1)
    closure = float(np.max(abs(solution.y[:d, -1] - state)))
    control = float(np.max(abs(matrix - shooting.y[d:, -1].reshape(d, d))))
    assert closure < 1e-8 and control < 1e-7
    item = dict(index=len(calls) + 1, rho=float(rho), period=period, residual=residual,
                multiplier=dict(real=float(dominant.real), imag=float(dominant.imag), modulus=float(abs(dominant))),
                closure=closure, strict_matrix_difference=control, shooting_history=history)
    filename = f'point-{len(calls)+1:02d}.npz'
    np.savez_compressed(out / filename, initial_state=state, parameters=p, weights=weights,
                        period=period, quotient_monodromy=matrix)
    item['npz_path'], item['npz_sha256'] = filename, sha(out / filename)
    calls.append(item)
    seeds.append((rho, state, period))
    write(out / 'analysis.json', record)
    print(json.dumps(item), flush=True)
    return residual


try:
    threshold = brentq(evaluate, below_rho['rho'], above_rho['rho'], xtol=5e-11, rtol=1e-12, maxiter=20)
    closest = min(calls, key=lambda p: abs(p['rho'] - threshold))
    lower, upper = threshold - 1e-5, threshold + 1e-5
    lower_value, upper_value = evaluate(lower), evaluate(upper)
    multiplier = complex(closest['multiplier']['real'], closest['multiplier']['imag'])
    symmetry = []
    permutation = [3, 0, 4, 1, 2]
    for point in (points[0], points[-1]):
        p = np.array([point['rho'], *branch['fixed_r_b']])
        x = np.array(point['initial_state'])
        solution = cycles.integrate(x, point['period'] / 6, p, weights, rtol=2e-13, max_step=.01)
        initial_rows = np.column_stack(cycles.unpack(x, weights))
        final_rows = np.column_stack(cycles.unpack(solution.y[:d if (d:=len(x)) else 0, -1], weights))
        error = float(np.max(abs(final_rows - initial_rows[permutation])))
        assert error < 1e-9
        symmetry.append(dict(rho=point['rho'], fraction='1/6', current_from_initial_group_permutation=permutation,
                             maximum_component_residual=error))
    record.update(status='complete', numerical_threshold=threshold, period_at_threshold=closest['period'],
                  multiplier_at_threshold=closest['multiplier'],
                  angle_radians=float(np.angle(multiplier)), angle_over_two_pi=float(np.angle(multiplier) / (2 * np.pi)),
                  low_order_resonance_distances={str(k): float(abs(multiplier ** k - 1)) for k in range(1, 5)},
                  centered_modulus_slope=(upper_value - lower_value) / (upper - lower),
                  slope_step=1e-5, numerical_sign_bracket=dict(lower=lower, lower_value=lower_value, upper=upper, upper_value=upper_value),
                  sixth_period_symmetry_controls=symmetry, wall_seconds=time.perf_counter() - started)
except Exception as error:
    record.update(status='failed', error=repr(error), wall_seconds=time.perf_counter() - started)
    raise
finally:
    write(out / 'analysis.json', record)

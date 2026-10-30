"""Short packaging controls for later motor modules in a clean public clone.

Uses an explicitly supplied clone and output file; no archived result is written.
Run with isolated Python from outside the source and destination repositories.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.integrate import solve_ivp

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
root = args.root.resolve()
assert not (root/'chua_gpu.py').exists()
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=root)
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
remote = subprocess.check_output(['git', 'remote', 'get-url', 'origin'], cwd=root, text=True).strip()
assert remote == 'https://github.com/Griboedow/chaos-attractor-lab.git'
sys.path.insert(0, str(root/'research'))
import motor_rows as model
import motor_rows_gpu as gpu_module
import motor_rows_activity_path as activity
import analyze_motor_cluster_cycle as cluster

records = []
for n in (3, 8):
    start, _ = model.source_start(n)
    parameters = np.array([2.934090909090909, 4.7551460477513405, 1.])
    tangent = gpu_module.default_tangent(n)
    gpu = gpu_module.MotorRowsGPU(n)
    result = gpu.scan(parameters[None], start[None], transient=0., duration=.2,
                      rtol=1e-10, step_factor=.5, max_wall=30., time_mode='original')
    d = model.dimension(n)
    def field(t, z):
        f, j = model.field_jac(z[:d], parameters, n)
        return np.r_[f, j@z[d:]]
    cpu = solve_ivp(field, (0, .2), np.r_[start, tangent], method='DOP853',
                    rtol=1e-12, atol=1e-14, max_step=.01)
    state_error = float(np.max(abs(result['final_states'][0]-cpu.y[:d, -1])))
    growth_error = float(abs(result['lle'][0]*.2-np.log(np.linalg.norm(cpu.y[d:, -1]))))
    assert cpu.success and result['status'][0] == 0
    assert max(state_error, growth_error) < 1e-8
    records.append(dict(rows=n, duration_original_tau=.2, device=gpu.device.name,
                        maximum_state_error=state_error, log_growth_error=growth_error,
                        passed=True))
assert 'chua_gpu' not in sys.modules
for module in (model, gpu_module, activity, cluster):
    assert Path(module.__file__).resolve().is_relative_to(root)
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=root)
record = dict(checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              verified_public_commit=revision, repository=remote, passed=True,
              isolated_python=bool(sys.flags.isolated), no_chua_import=True,
              fresh_public_clone_environment_reused_from_initial_validation=True,
              new_install_performed_for_this_check=False, tracked_evidence_unchanged=True,
              controls=records, control_driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              scope='Short N3 and N8 GPU/CPU packaging controls and activity/cluster module imports after updating the independent public clone. Not a new science survey or validation of every archived trajectory.')
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
print(json.dumps(record, indent=2))

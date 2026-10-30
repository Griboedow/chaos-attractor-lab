"""Independent public-clone packaging checks; writes only to an explicit output."""
import argparse
import contextlib
import datetime
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.integrate import solve_ivp

parser = argparse.ArgumentParser()
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
root = args.root.resolve()
out = args.output.resolve()
assert (root / '.git').is_dir()
assert not (root / 'chua_gpu.py').exists()
out.mkdir(parents=True, exist_ok=True)
before = subprocess.check_output(['git', 'status', '--porcelain'], cwd=root)
assert not before, before
commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
remote = subprocess.check_output(['git', 'remote', 'get-url', 'origin'], cwd=root, text=True).strip()
assert remote == 'https://github.com/Griboedow/chaos-attractor-lab.git'
sys.path.insert(0, str(root / 'research'))
records = {}
with (out / 'runtime.log').open('w', encoding='utf8') as log, contextlib.redirect_stdout(log):
    import epithelial_chain as ec
    import rheostat as rh
    import photofeedback as pf
    from thermal_sheet_gpu import ThermalSheetGPU
    from thermomechanical_sheet import field_jacobian as sheet_field
    chain = ec.ChainGPU(4)
    ec.validate(chain, out / 'chain')
    rheostat = rh.RheostatGPU()
    rh.validate(rheostat, out / 'rheostat')
    photo = pf.PhotoFeedbackGPU()
    thermal = ThermalSheetGPU()
    for name, gpu in [('chain', chain), ('rheostat', rheostat), ('photo', photo), ('thermal', thermal)]:
        records[name] = {'device': gpu.info['name'], 'precision': gpu.info['precision'], 'constructor_passed': True}
    p, s = pf.parameter_at_equilibrium(1., .1, 1., 1., .2)
    s = s + np.array([.01, 0, 0])
    specs = [('photo', photo, p, s, lambda z, p: pf.rhs(0, z, p)),
             ('thermal', thermal, np.array([1., 1., 1., .1, .1]), np.array([.5, .01, .5]),
              lambda z, p: np.r_[sheet_field(z[:3], p)[0], sheet_field(z[:3], p)[1] @ z[3:]])]
    for name, gpu, p, s, rhs in specs:
        result = gpu.scan(p[None], s[None], transient=0, duration=1, rtol=1e-10, atol=1e-13)
        cpu = solve_ivp(lambda t, z: rhs(z, p), (0, 1), np.r_[s, np.array([1., 2., 3.]) / np.sqrt(14)],
                        method='DOP853', rtol=1e-12, atol=1e-14)
        state_error = float(max(abs(result['final_states'][0] - cpu.y[:3, -1])))
        lle_error = float(abs(result['lle'][0] - np.log(np.linalg.norm(cpu.y[3:, -1]))))
        passed = bool(cpu.success and result['status'][0] == 0 and max(state_error, lle_error) < 1e-7)
        records[name].update(state_error=state_error, lle_error=lle_error, control_passed=passed)
        assert passed, (name, state_error, lle_error)
    for name in ['chain', 'rheostat']:
        records[name]['validation'] = json.loads((out / name / 'validation.json').read_text())
    import dcr
    from sksundae.ida import IDA
    ref = np.load(root / 'results/dcr/ida-analytic-1-fine.npz')
    p = dcr.Parameters()
    y0 = ref['initial_state']
    def residual(t, y, yp, residual_out):
        residual_out[:] = dcr.residual(y, yp, p)
    def jac(t, y, yp, r, cj, jac_out):
        jac_out[:] = dcr.residual_jacobian(y, yp, cj, p)
    sol = IDA(residual, jacfn=jac, rtol=1e-9, atol=1e-9, max_step=2.5e-5).solve(
        np.linspace(0, .01, 11), y0, dcr.field(y0, p))
    err = float(np.max(abs(sol.y - ref['states'][:11])))
    assert sol.success and err < 1e-7
    records['dcr'] = {'duration': .01, 'maximum_state_error_against_archived_fine': err,
                      'passed': True, 'scope': 'Short optional-dependency packaging control only'}
    assert 'chua_gpu' not in sys.modules
    for module in [ec, rh, pf, dcr]:
        assert Path(module.__file__).resolve().is_relative_to(root)

pip = subprocess.run([sys.executable, '-m', 'pip', 'check'], capture_output=True, text=True, check=True)
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=root)
result = dict(passed=True, checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              repository=remote, verified_commit=commit, isolated_python=bool(sys.flags.isolated),
              fresh_public_clone=True, fresh_venv=True, pip_check=pip.stdout.strip(),
              dependencies={n: importlib.metadata.version(n) for n in ['numpy','scipy','matplotlib','pyopencl','pillow','scikit-sundae']},
              controls=records, no_chua_module_imported=True, original_evidence_modified=False,
              control_driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              scope='Four extracted-runtime GPU/CPU controls and one short optional DCR control in a fresh public clone and new environment, executed from outside both repositories. This is packaging verification, not a new parameter survey or validation of every archived model.')
(out / 'public-clone-validation.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf8')
print(json.dumps(result, indent=2))

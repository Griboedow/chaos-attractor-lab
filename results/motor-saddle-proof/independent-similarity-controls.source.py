"""Independent similarity controls with known real and complex roots."""
from __future__ import annotations
import datetime
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'research'))
import flint
from flint import arb, arb_mat, ctx
from certify_motor_saddle_spectrum import similarity_schur
from validate_motor_cycle import ball_record, read_ball

OUTPUT = ROOT/'results/motor-saddle-proof/independent-similarity-controls'
if OUTPUT.exists():
    raise FileExistsError('Preserve earlier controls')
OUTPUT.mkdir()
ctx.prec = 256

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def basis_file(name, entries):
    path = OUTPUT/name
    path.write_text(json.dumps(dict(basis=[[[str(v.numerator), str(v.denominator)]
        for v in row] for row in entries]), indent=2)+'\n', encoding='utf-8')
    return path

entries = [[Fraction(int(i == j)) for j in range(5)] for i in range(5)]
entries[1][1], entries[2][2] = Fraction(1, 256), Fraction(256)
for i in range(4):
    entries[i][i+1] = Fraction(64)
path = basis_file('nonorthogonal-exact-basis.json', entries)
q = arb_mat([[arb(v.numerator)/v.denominator for v in row] for row in entries])
qi = q.inv()
assert q.det() == arb(1)
assert not q.transpose() == qi
checks, cases = [], []

def target(scalar):
    # Eigenvalues: 1 (known flow factor), scalar, -1/4, 1/4 +/- i/2.
    a = arb_mat(5, 5)
    a[0,0] = 1
    a[1,1] = scalar
    a[2,2] = arb('-1/4')
    a[3,3] = a[4,4] = arb('1/4')
    a[3,4], a[4,3] = arb('-1/2'), arb('1/2')
    return a

for name, scalar, expected in (
    ('stable_real_and_complex_roots', arb('1/2'), True),
    ('stable_small_interval_uncertainty', arb('1/2', '1e-35'), True),
    ('wide_uncertainty_includes_expanding_root', arb('1/2', '0.6'), False),
    ('exact_expanding_root', arb('5/4'), False),
    ('extra_exact_unit_root', arb(1), False)):
    a = target(scalar)
    matrix = q*a*qi
    result = similarity_schur(matrix, path, arb('0.8'))
    transformed = arb_mat([[read_ball(v) for v in row] for row in result['similar_monodromy']])
    enclosure = all(transformed[i,j].contains(a[i,j]) for i in range(5) for j in range(5))
    retained = scalar.rad().is_zero() or not transformed[1,1].rad().is_zero()
    check = dict(name=name, passed=bool(result['passed'] == expected and enclosure and retained),
        expected_spectral_pass=expected, actual_spectral_pass=result['passed'],
        target_enclosed=enclosure, interval_uncertainty_retained=bool(retained),
        target_scalar=ball_record(scalar), transformed_scalar=ball_record(transformed[1,1]),
        input_max_radius=ball_record(max(matrix[i,j].rad() for i in range(5) for j in range(5))))
    checks.append(check)
    cases.append(dict(name=name, certificate=result))

singular = [list(row) for row in entries]
singular[-1] = list(singular[0])
invalid = [list(row) for row in entries]
invalid[0][0] = Fraction(1, 3)
for name, entries2, exception in (
    ('singular_basis_rejected', singular, ArithmeticError),
    ('nondyadic_basis_rejected', invalid, ValueError)):
    rejected = False
    try:
        similarity_schur(target(arb('1/2')), basis_file(name+'.json', entries2), arb('0.8'))
    except exception:
        rejected = True
    checks.append(dict(name=name, passed=rejected))

# Verify the actual postprocessor dependency provenance and the source
# snapshots copied from the archived existence certificate.
actual_path = ROOT/'results/motor-saddle-proof/similarity-certificate.json'
actual = json.loads(actual_path.read_text())
dependency_matches = {name: sha(ROOT/'research'/name) == actual['source_sha256'][name]
    for name in ('validate_motor_cycle.py', 'validate_motor_paired_cycle.py')}
snapshot_matches = {name: sha(actual_path.with_name(actual_path.stem+'.'+name)) == digest
    for name, digest in actual['source_sha256'].items()}
checks.append(dict(name='actual postprocessor dependencies and archived source snapshots',
    dependency_matches=dependency_matches, snapshot_matches=snapshot_matches,
    passed=all(dependency_matches.values()) and all(snapshot_matches.values())))
record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    status='passed' if all(v['passed'] for v in checks) else 'failed',
    python=sys.version, python_flint=flint.__version__, precision_bits=ctx.prec,
    basis_determinant=ball_record(q.det()),
    basis_max_entry=ball_record(max(q[i,j].abs_upper() for i in range(5) for j in range(5))),
    inverse_max_entry=ball_record(max(qi[i,j].abs_upper() for i in range(5) for j in range(5))),
    basis_is_nonorthogonal=True, checks=checks, cases=cases,
    source_sha256={name:sha(ROOT/'research'/name) for name in
        ('certify_motor_saddle_spectrum.py','validate_motor_paired_cycle.py','validate_motor_cycle.py')},
    helper_sha256=sha(Path(__file__)), actual_certificate_sha256=sha(actual_path),
    limitation='Known-matrix implementation controls; these do not establish the motor orbit or replace its validated monodromy enclosure.')
(OUTPUT/'analysis.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
print(record['status'], len(checks), 'checks')
if record['status'] != 'passed':
    raise SystemExit(2)

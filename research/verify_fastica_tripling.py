"""Verify the retained structured FastICA packet without optional arrays."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json

ROOT=Path(__file__).resolve().parents[1];PACKET=ROOT/'results/fastica-tripling'
def read(name):return json.loads((PACKET/name).read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
manifest=read('manifest.json')
for row in manifest['files']:
    path=PACKET/row['path'];assert path.stat().st_size==row['bytes'] and sha(path)==row['sha256'],row['path']
model=sha(PACKET/'code/model.py')
for protocol in ['pilot-01/protocol.json','long-01/protocol.json']:
    row=read(protocol)
    if protocol.startswith('pilot'):assert row['source_sha256']==model
    else:assert row['model_sha256']==model and row['source_sha256']==sha(PACKET/'code/long_check.py')
e=Q(1,25);k=Q(2,5);L=1-18*e*k-36*e**2*k**2-18*e**3*k**3;U=6*e+36*e**2*k+54*e**3*k**2+108*e**4*k**3
exact=read('exact-check.json');assert exact['status']=='passed'
assert L==Q(exact['L']) and U==Q(exact['U']) and Q(exact['ratio_bound'])==U/L<k and Q(exact['strict_margin'])==k*L-U>0
assert exact['source_sha256']==sha(PACKET/'code/exact_check.py') and exact['model_sha256']==model
pilot=read('pilot-01/report.json');assert len(pilot['rows'])==80
assert pilot['positive_cases']==list(range(1,36))+[41]
assert all(r['status']==0 and r['completed']==105000 for r in pilot['rows'])
long=read('long-01/report.json');assert len(long['rows'])==4
assert all(r['status']==0 and r['completed']==5010000 and len(r['rates'])==10 and min(x[0] for x in r['rates'])>0 for r in long['rows'])
fresh=read('long-01/fresh-starts.json');assert len(fresh)==2
for row in fresh:
    assert len(row['starts'])==128 and row['positive_cases']==list(range(128))
    assert all(st==[0,210000] for st in row['status'])
    assert all(len(rates)==4 and min(r[0] for r in rates)>.005 for rates in row['rates'])
precision=read('precision.json');assert len(precision)==4
for row in precision:
    assert row['source_sha256']==sha(PACKET/'code/precision_check.py')
    assert row['observed']==20000 and len(row['rates'])==5 and min(map(float,row['rates']))>0
library=read('library.json');assert len(library)==2
for row in library:
    assert row['source_sha256']==sha(PACKET/'code/library_check.py') and row['model_sha256']==model
    assert row['count']==52000 and row['window_counts']==[10000]*5 and not row['first_unit_finished']
    assert min(row['rates'])>0 and row['max_polynomial_step_error']<4e-13
controls=read('controls.json');assert controls['source_sha256']==sha(PACKET/'code/controls.py')
assert len(controls['empirical_controls'])==8
assert max(r['maximum_step_difference'] for r in controls['empirical_controls'])<1e-11
assert max(r['maximum_tangent_relative_difference'] for r in controls['empirical_controls'])<7e-12
omitted=read('omitted-data.json');registry=json.loads((ROOT/'records/generated-data.json').read_text())
table={r['path']:r for r in registry['files']}
for row in omitted:assert table[row['path']]==row
print(json.dumps({'status':'passed','packet_files':len(manifest['files']),'pilot_cases':80,'long_runs':4,'fresh_starts':256,'precision_runs':4,'library_runs':2,'omitted_bytes':sum(r['bytes'] for r in omitted),'topological_proof':'Written analytic argument; not machine formalized.'}))

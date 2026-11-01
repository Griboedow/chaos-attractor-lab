"""Audit frozen evidence, selection, inputs, controls and unthinned counts."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2];P=ROOT/'results/fastica'
read=lambda p:json.loads((P/p).read_text(encoding='utf-8'))
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=read('manifest.json')
for row in manifest['files']:
    p=ROOT/row['path'];assert p.stat().st_size==row['bytes'] and digest(p)==row['sha256'],row['path']
pilot=read('pilot-01/report.json');assert len(pilot['cases'])==576
retained=[r['case'] for r in pilot['cases'] if min(r['rates'])>.005 and r['period_up_to_sign']==0 and r['minimum_update_norm']>1e-8 and r['status']==0]
assert retained==pilot['retained']==[15,40,204,246,351]
for case in retained:
    d=read(f'follow-01/case-{case}-input.json');X=np.array(d['X']);s=pilot['cases'][case]
    assert X.shape==(s['n'],4)
    assert hashlib.sha256(X.tobytes()).hexdigest()==s['input_sha256']
    assert np.max(abs(X.mean(0)))<1e-13 and np.max(abs(X.T@X/len(X)-np.eye(4)))<1e-12
long=read('long-01/report.json');assert [r['case'] for r in long]==retained
for r in long:
    a=np.array(r['rates']);assert a.shape==(10,3) and np.isfinite(a).all()
    assert np.allclose(a.mean(0),r['mean'],atol=1e-14,rtol=0)
    assert r['observed']==5000000 and r['minimum_update_norm']>1e-8 and r['minimum_qr_diagonal']>0
    if r['case']!=40:assert min(a[:,0])>.005
    else:assert np.max(a[-1])<0
basin=read('basin-01/report.json')['rows'];assert len(basin)==320
counts={str(c):sum(r['retained'] for r in basin if r['case']==c) for c in retained}
assert counts=={'15':10,'40':0,'204':64,'246':64,'351':10}
perturb=read('perturb-01/report.json');assert len(perturb)==96
assert [sum(r['retained'] for r in perturb if r['epsilon']==e) for e in [1e-6,1e-3,1e-2]]==[32,32,31]
for digits,cases in [(70,[15,204,246,351]),(100,[204,351])]:
    rows=read(f'precision-{digits}.json');assert [r['case'] for r in rows]==cases
    for r in rows:
        assert min(map(float,r['rates']))>.005 and len(r['rates'])==4
        assert r['input_sha256']==digest(P/f"follow-01/case-{r['case']}-input.json")
        assert r['source_sha256']==digest(P/'source-original/precision.py')
library=read('library-controls-v2.json');assert [r['case'] for r in library]==retained
for r in library:
    assert r['maximum_step_identity_error']<2e-13
    assert r['source_probe_sha256']==digest(P/'source-original/library_probe.py')
    for n,v in zip(r['window_counts'],r['complete_window_rates']):assert (v is not None)==(n==10000)
    if r['case'] in [15,204,351]:assert min(r['complete_window_rates'])>.005 and not r['first_unit_converged']
controls=read('controls.json');assert controls['maximum_step_error']<1e-12 and controls['maximum_tangent_relative_error']<2e-7
assert max(controls['outcomes'][0]['return_multiplier_moduli'])<1
assert min(controls['outcomes'][1]['continuation_rates'])>.005
stopping=read('stopping.json');r=next(r for r in stopping['rows'] if r['case']==204)
assert r['iterations_before_first_unit_exit']==6 and min(r['continuation_rates'])>.005
for r in read('histograms/render-inputs.json')['rows']:
    p=P/f"histograms/case-{r['case']}-histogram.npz";assert digest(p)==r['histogram_sha256']
    with np.load(p) as data:assert int(data['histogram'].sum())==r['state_count']==1000000
    assert r['clipped']==0 and r['columns_zero_based']==[0,1,3]
for r in read('portability.json'):
    assert digest(P/'source-original'/r['file'])==r['original_sha256']
    assert digest(ROOT/'research/fastica'/r['file'])==r['portable_sha256']
print(json.dumps({'manifest_files':len(manifest['files']),'screen_cases':576,'long_cases':5,'fresh_starts':320,'perturbations':96,'portrait_counts':[1000000,1000000],'checks':'passed'}))

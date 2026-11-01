"""Verify retained refinement evidence and optionally replay portable algebra."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
PACKET=ROOT/'results/hydrogel-variational-checks'

def read(path):return json.loads(path.read_text(encoding='utf-8'))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--numerical',action='store_true')
    ap.add_argument('--workdir',type=Path)
    args=ap.parse_args()
    manifest=read(PACKET/'manifest.json')
    for row in manifest['files']:
        p=PACKET/row['path']
        assert p.stat().st_size==row['bytes'] and digest(p)==row['sha256'],p
    for row in manifest['dependencies']:assert digest(ROOT/row['path'])==row['sha256'],row
    runs={p.stem:read(p) for p in (PACKET/'runs').glob('*.json')}
    workflow=read(PACKET/'workflow.json')
    assert len(runs)==20 and {r['run'] for r in workflow}==set(runs)
    assert sum(r['success'] for r in runs.values())==19
    assert not runs['field130-continuation']['success']
    arrays=read(PACKET/'arrays.json')
    assert len(arrays['files'])==40
    hashes={a['sha256'] for a in arrays['files']}
    hashes.update(a['trajectory']['sha256'] for a in arrays['initial_inputs'])
    sources={digest(p) for p in (PACKET/'code').glob('*.py')}
    for name,r in runs.items():
        assert r['code_sha256'] in sources and r['solver_sha256'] in sources,name
        key='seed_checkpoint' if 'seed_checkpoint' in r else 'seed'
        assert r[key]['sha256'] in hashes,name
        assert r['completed_time']<=r['final_accepted_time']+1e-12
        if r['success']:
            assert abs(r['final_accepted_time']-r['time_origin']-r['requested_duration'])<1e-9
        if 'observed_duration' in r:
            initial=r.get('initial_observed_duration',0)
            expected=initial+r['final_accepted_time']-r['time_origin']-r['burn_duration']
            assert abs(r['observed_duration']-expected)<1e-9,name
            assert len(r['exponents'])==3 and r['observed_duration']>0
    extended=runs['spectrum124-extension']
    assert extended['initial_observed_duration']==runs['spectrum124-n64']['observed_duration']==96
    assert extended['observed_duration']==336
    clearance=read(PACKET/'analysis/clearance-later.json')
    failure=next(r for r in clearance if r['run']=='field130-continuation')
    assert failure['samples_below_threshold']==496 and failure['samples']==2029
    assert failure['minimum_clearance']<failure['tube_diameter_threshold']
    drift=read(PACKET/'analysis/drift.json')
    windows=next(r for r in drift['runs'] if r['run']=='spectrum124-extension')['windows']
    assert windows[0]['mean'][2]-windows[-1]['mean'][2]>.03
    checked=[]
    if args.numerical:
        if args.workdir is None:ap.error('--numerical requires a new external --workdir')
        target=args.workdir.resolve()
        if target.is_relative_to(ROOT):ap.error('Keep numerical outputs outside the checkout')
        target.mkdir(parents=True,exist_ok=False)
        for p in (PACKET/'code').glob('*.py'):shutil.copyfile(p,target/p.name)
        env=os.environ.copy()
        env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONUTF8='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_CACHE_DIR=str(target/'numba'),TEMP=str(target),TMP=str(target))
        commands=[['check_portable.py'],['check_tangent.py'],['check_history.py'],
                  ['recurrence.py','--check','--output',str(target/'recurrence-check.json')]]
        for command in commands:
            result=subprocess.run([sys.executable,*command],cwd=target,env=env,check=True,capture_output=True,text=True)
            (target/(Path(command[0]).stem+'-stdout.txt')).write_text(result.stdout,encoding='utf-8')
            checked.append(command[0])
    print(json.dumps(dict(retained_hashes=len(manifest['files'])+len(manifest['dependencies']),
        runs=len(runs),complete=19,failed=1,external_arrays=len(arrays['files']),numerical_checks=checked,
        scope='Retained hashes, run/seed accounting and optional synthetic algebra. No long integration replay, continuum convergence or chaotic-attractor certificate.')))

if __name__=='__main__':main()

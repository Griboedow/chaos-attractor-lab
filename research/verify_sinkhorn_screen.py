"""Check the retained Sinkhorn campaign and optionally replay small-model checks."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PACKET=ROOT/'results/sinkhorn-screen'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--numerical',action='store_true');ap.add_argument('--workdir',type=Path)
    args=ap.parse_args();manifest=read(PACKET/'manifest.json')
    for r in manifest['files']:
        p=PACKET/r['path'];assert p.stat().st_size==r['bytes'] and digest(p)==r['sha256'],p
    for r in manifest['dependencies']:assert digest(ROOT/r['path'])==r['sha256'],r
    code={p.name:digest(p) for p in (PACKET/'code').glob('*.py')}
    originals={'/'.join(r['path'].replace('\\','/').split('/')[-2:]):r for r in read(PACKET/'original-summaries.json')['files']}
    flagged=set();all_cases=[]
    for n in (2,3,4):
        pilot=read(PACKET/'pilot'/f'n{n}.json');assert len(pilot['runs'])==192
        assert pilot['source_sha256']=={name:code[name] for name in ('survey.py','model.py')}
        for r in pilot['runs']:
            all_cases.append(r);assert np.array(r['log_kernel']).shape==(n,n)
            assert len(r['initial_state'])==2*n-2 and 0<r['omega']<2
            if not r['capture_step'] and r['finite_rate']>1e-5:flagged.add((n,r['case']))
    assert len(all_cases)==576 and sum(bool(r['capture_step']) for r in all_cases)==314 and len(flagged)==40
    for stage,start,block in [('followup',30000,25000),('extension',530000,250000)]:
        actual=set()
        for n in (3,4):
            idx=read(PACKET/stage/f'n{n}-index.json');assert idx['block_length']==block and idx['blocks']==20
            assert all(code[k]==v for k,v in idx['sources'].items())
            if stage=='extension':assert idx['parent_sha256']==originals[f'followup01/n{n}.json']['sha256']
            for case in idx['cases']:
                r=read(PACKET/stage/f'n{n}-case{case}.json');actual.add((n,case))
                assert r['start_step']==start and len(r['blocks'])==20
                assert [b['stop_step'] for b in r['blocks']]==[start+block*(i+1) for i in range(20)]
                if stage=='extension':
                    previous=read(PACKET/'followup'/f'n{n}-case{case}.json')
                    assert r['starting_state']==previous['blocks'][-1]['checkpoint']
        assert actual==flagged
    recurrence=read(PACKET/'analysis/recurrence.json');assert recurrence['source_sha256']==code['analyze.py']
    rows=recurrence['runs'];assert {(r['shape'][0],r['case']) for r in rows}==flagged
    assert sum(r['last_block_rate']>.01 for r in rows)==11
    assert sum(r['smallest_tested_recurrence']==1 for r in rows)==4
    assert sum(r['smallest_tested_recurrence'] is not None and r['smallest_tested_recurrence']>1 for r in rows)==21
    spectra=read(PACKET/'analysis/spectra.json');assert len(spectra['runs'])==4
    for r in spectra['runs']:
        rates=np.array(r['rates']);assert rates.shape==(20,sum(r['shape'])-2)
        assert np.max(abs(rates.sum(axis=1)-r['expected_volume_rate']))<1e-10
        assert np.max(abs(rates.mean(axis=0)-r['mean_rates']))<1e-14
    controls=read(PACKET/'analysis/controls.json');assert controls['source_sha256']==code['controls.py'] and len(controls['runs'])==24
    assert sum(r['rates'] is None for r in controls['runs'])==4
    for folder,parent in [('precision01','followup01/n4.json'),('precision02','extension01/n4.json')]:
        meta=read(PACKET/folder/'metadata.json');assert meta['source_sha256']==code['precision.py']
        assert meta['parent_sha256']==originals[parent]['sha256']
        runs=[read(p) for p in (PACKET/folder).glob('*digits.json')]
        assert {r['precision_decimal_digits'] for r in runs}=={70,100}
        assert all(r['observed']==20000 and min(r['rates'])>0 for r in runs)
    arrays=read(PACKET/'arrays.json');assert len(arrays['files'])==11
    assert sum(r['bytes'] for r in arrays['files'])==arrays['total_bytes']==46513332
    replay=[]
    if args.numerical:
        if args.workdir is None:ap.error('--numerical requires a new external --workdir')
        work=args.workdir.resolve()
        if work.is_relative_to(ROOT):ap.error('Keep replay outputs outside the checkout')
        work.mkdir(parents=True,exist_ok=False)
        for name in ('model.py','check_model.py'):shutil.copyfile(PACKET/'code'/name,work/name)
        env=os.environ.copy();env.update(PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',NUMBA_CACHE_DIR=str(work/'numba'),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',TEMP=str(work),TMP=str(work))
        r=subprocess.run([sys.executable,'check_model.py'],cwd=work,env=env,check=True,capture_output=True,text=True)
        (work/'stdout.txt').write_text(r.stdout,encoding='utf-8');replay.append('18 independent small-model checks')
    print(json.dumps(dict(retained_hashes=len(manifest['files'])+len(manifest['dependencies']),pilot_cases=576,flags=40,continuation_endpoint=5530000,
        final_positive=11,nontrivial_recurrences=21,fixed_state_recurrences=4,other_weak_or_drifting=4,external_bytes=arrays['total_bytes'],numerical_replay=replay,
        scope='Retained evidence and accounting; no infinite-time chaos, basin, or historical-novelty certificate.')))

if __name__=='__main__':main()

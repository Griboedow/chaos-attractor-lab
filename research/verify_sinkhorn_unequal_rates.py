"""Verify compact unequal-rate evidence and optionally replay independent checks."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1];PACKET=ROOT/'results/sinkhorn-unequal-rates'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--numerical',action='store_true');ap.add_argument('--workdir',type=Path);args=ap.parse_args()
    manifest=read(PACKET/'manifest.json')
    for r in manifest['files']:
        p=PACKET/r['path'];assert p.stat().st_size==r['bytes'] and digest(p)==r['sha256'],p
    code={p.name:digest(p) for p in (PACKET/'code').glob('*.py')};datasets={};originals={}
    for entry in read(PACKET/'datasets.json')['datasets']:
        index=read(PACKET/entry['index']);rows=[r for part in index['parts'] for r in read(PACKET/part)['runs']];assert len(rows)==index['runs'];d=index['metadata']|{'runs':rows};datasets[entry['original']]=d;originals[entry['original']]=entry['original_sha256'];assert index['original_sha256']==entry['original_sha256']
        for k,v in d['sources'].items():assert code[k]==v,(entry,k)
    assert len(datasets)==26
    for stage,shapes,count,expected,origin in [('pilot01',[(2,3),(3,2),(2,4),(4,2)],384,20,60000),('dual01',[(2,2),(2,3),(3,2)],768,12,30000),('balanced01',[(3,2)],1024,57,40000)]:
        flags=set()
        for m,n in shapes:
            rows=datasets[f'{stage}/{m}x{n}.json']['runs'];assert len(rows)==count and [r['case'] for r in rows]==list(range(count))
            for r in rows:
                assert np.array(r['log_kernel']).shape==(m,n) and len(r['initial_state'])==m+n-2
                flagged=(not r['capture_step'] and r['finite_rate']>1e-5) if stage=='pilot01' else bool(r['rates']) and r['rates']['mean']>.005
                if flagged:flags.add((m,n,r['case']))
                if stage!='pilot01':assert abs(r['determinant']-(1-r['wx'])**(m-1)*(1-r['wy'])**(n-1))<1e-14
        assert len(flags)==expected
        first={'pilot01':'followup01','dual01':'dual-follow01','balanced01':'balanced-follow01'}[stage]
        actual={(m,n,r['case']) for m,n in shapes for r in datasets[f'{first}/{m}x{n}.json']['runs']};assert actual==flags
        for m,n in shapes:
            child=datasets[f'{first}/{m}x{n}.json'];assert child['parent_sha256']==originals[f'{stage}/{m}x{n}.json']
            for r in child['runs']:
                assert all(b['stop_step']==origin+(i+1)*child['block_length'] for i,b in enumerate(r['blocks']))
    for parent,child,shapes in [('followup01','rectangular-long01',[(2,3),(3,2),(2,4),(4,2)]),('dual-follow01','dual-long01',[(2,2),(2,3),(3,2)]),('balanced-follow01','balanced-long01',[(3,2)]),('balanced-long01','selected-100m',[(3,2)])]:
        for m,n in shapes:
            a=datasets[f'{parent}/{m}x{n}.json'];b=datasets[f'{child}/{m}x{n}.json'];assert b['parent_sha256']==originals[f'{parent}/{m}x{n}.json']
            eligible={r['case'] for r in a['runs'] if r.get('completed_steps',a['block_length']*a['blocks'])==a['block_length']*a['blocks']}
            if child=='selected-100m':eligible&={150,987}
            assert {r['case'] for r in b['runs']}==eligible
            for r in b['runs']:
                previous=next(v for v in a['runs'] if v['case']==r['case']);assert r['starting_state']==previous['blocks'][-1]['checkpoint']
    final=datasets['selected-100m/3x2.json'];assert final['stop_step']==106040000
    assert all(r['completed_steps']==100000000 and len(r['blocks'])==20 and min(b['rate'] for b in r['blocks'])>0 for r in final['runs'])
    analysis=read(PACKET/'analysis.json');assert analysis['source_sha256']==code['analyze_runs.py']
    assert [r['inputs'] for r in analysis['campaigns']]==[1536,2304,1024]
    assert [r['positive_flags'] for r in analysis['campaigns']]==[20,12,57]
    assert [r['noncontracting_volume'] for r in analysis['campaigns']]==[0,76,0]
    balanced=[r for r in analysis['final_windows'] if r['stage']=='balanced-long01'];assert sum(r.get('last_growth',-1)>0 for r in balanced)==51
    for filename in ('spectra-selected-106m.json','spectra-balanced-1m.json','spectra-rectangular-6m.json'):
        d=read(PACKET/filename);assert d['parent_sha256'] in originals.values()
        for k,v in d['sources'].items():assert code[k]==v
        for r in d['runs']:
            assert r['completed'];a=np.array(r['rates']);assert a.shape==(20,sum(r['shape'])-2)
            expected=(r['shape'][0]-1)*np.log(abs(1-r['wx']))+(r['shape'][1]-1)*np.log(abs(1-r['wy']))
            assert abs(expected-r['expected_volume_rate'])<1e-14 and max(abs(a.sum(1)-expected))<1e-10
            assert max(abs(a.mean(0)-r['mean_rates']))<1e-14
            if filename=='spectra-selected-106m.json':assert sum(a.mean(0)>0)==2
    for folder,observed in [('precision-balanced01',20000),('precision-rectangular01',20000),('precision-balanced-6m',200000)]:
        meta=read(PACKET/folder/'metadata.json');assert meta['source_sha256']==code['precision_dual.py'] and meta['parent_sha256'] in originals.values()
        for p in (PACKET/folder).glob('*digits.json'):
            r=read(p);assert r['completed'] and r['observed']==observed and len(r['rates'])*2000==observed and min(r['rates'])>0
    for filename,key in [('schur-check.json','check_schur.py'),('actual-checks.json','check_actual_dual.py')]:
        d=read(PACKET/filename)
        if 'source_sha256' in d:assert d['source_sha256']==code[key]
        else:
            for k,v in d['sources'].items():assert code[k]==v
    assert len(read(PACKET/'schur-check.json')['checks'])==240 and len(read(PACKET/'actual-checks.json')['checks'])==18
    arrays=read(PACKET/'arrays.json');assert sum(r['bytes'] for r in arrays['files'])==arrays['total_bytes']
    replay=[]
    if args.numerical:
        if args.workdir is None:ap.error('--numerical requires --workdir')
        work=args.workdir.resolve()
        if work.is_relative_to(ROOT):ap.error('Keep numerical outputs outside the checkout')
        work.mkdir(parents=True,exist_ok=False)
        for p in (PACKET/'code').glob('*.py'):shutil.copyfile(p,work/p.name)
        for name,d in datasets.items():
            p=work/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8')
        env=os.environ.copy();env.update(PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',NUMBA_CACHE_DIR=str(work/'numba'),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',TEMP=str(work),TMP=str(work))
        for name in ('check_schur.py','check_actual_dual.py'):
            r=subprocess.run([sys.executable,name],cwd=work,env=env,capture_output=True,text=True,check=True);(work/(name+'.stdout.txt')).write_text(r.stdout,encoding='utf-8');replay.append(name)
    print(json.dumps(dict(retained_files=len(manifest['files']),pilot_inputs=4864,rectangular_flags=20,unequal_flags=12,balanced_flags=57,selected_endpoint=106040000,selected_positive_qr_rates=2,external_array_bytes=arrays['total_bytes'],numerical_replay=replay,scope='Numerical evidence and structural checks; no asymptotic attraction or historical-novelty certificate.')))

if __name__=='__main__':main()

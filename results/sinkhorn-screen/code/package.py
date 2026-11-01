"""Retain inputs and summaries; describe every externally stored state array."""
from pathlib import Path
import argparse,hashlib,json,platform,shutil
import numpy as np,scipy,numba,mpmath,matplotlib

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,value,compact=False):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(value,indent=None if compact else 2,separators=(',',':') if compact else None,allow_nan=False)+'\n',encoding='utf-8')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo',type=Path,required=True);args=ap.parse_args()
    target=args.repo/'results/sinkhorn-screen';target.mkdir(parents=True,exist_ok=False)
    code=['model.py','survey.py','followup.py','extend.py','precision.py','spectrum.py','analyze.py','controls.py','check_model.py','check_actual.py','diagnostic.py','preview.py','cache_sources.py','package.py']
    (target/'code').mkdir()
    for name in code:shutil.copyfile(name,target/'code'/name)
    originals=[]
    for stage in ('pilot02','followup01','extension01'):
        for source in sorted(Path(stage).glob('n*.json')):
            data=json.loads(source.read_text());n=data.get('size',data['runs'][0]['shape'][0]);folder=stage.replace('02','').replace('01','')
            if stage=='pilot02':
                dst=target/folder/source.name;dump(dst,data,compact=True);assert dst.stat().st_size<256*1024
                assert json.loads(dst.read_text())==data
                paths=[dst.relative_to(target).as_posix()]
            else:
                metadata={k:v for k,v in data.items() if k!='runs'};metadata['cases']=[r['case'] for r in data['runs']]
                dst=target/folder/f'n{n}-index.json';dump(dst,metadata);paths=[dst.relative_to(target).as_posix()]
                for r in data['runs']:
                    dst=target/folder/f"n{n}-case{r['case']}.json";dump(dst,r);paths.append(dst.relative_to(target).as_posix())
                assert [{**json.loads((target/p).read_text())} for p in paths[1:]]==data['runs']
            originals.append(dict(path=str(source.resolve()),sha256=digest(source),bytes=source.stat().st_size,retained=paths,conversion='JSON values preserved; whitespace compacted or runs split into separate files.'))
    for source,dest in [('model-check.json','checks/model.json'),('actual-check.json','checks/actual-states.json'),
                        ('recurrence02.json','analysis/recurrence.json'),('controls02.json','analysis/controls.json'),
                        ('spectra02/summary.json','analysis/spectra.json'),('sources.json','sources.json')]:
        dst=target/dest;dst.parent.mkdir(exist_ok=True,parents=True);shutil.copyfile(source,dst)
    for folder in ('precision01','precision02'):
        for source in Path(folder).glob('*.json'):
            dst=target/folder/source.name;dst.parent.mkdir(exist_ok=True,parents=True);shutil.copyfile(source,dst)
    arrays=[]
    for folder in ('pilot02','followup01','extension01','precision01','precision02'):
        for source in sorted(Path(folder).glob('*')):
            if source.suffix not in ('.npy','.npz'):continue
            data=np.load(source,allow_pickle=False)
            shapes={k:list(data[k].shape) for k in data.files} if source.suffix=='.npz' else {'states':list(data.shape)}
            arrays.append(dict(path=str(source.resolve()),bytes=source.stat().st_size,sha256=digest(source),shapes=shapes,
                coverage='Endpoints and the final 512 pilot steps; early-captured tails contain NaN outside their declared coverage.' if folder=='pilot02' else 'Every state of the final 25,000 steps of this continuation, not the full preceding trajectory.' if folder in ('followup01','extension01') else 'Every observed state of the 20,000-step precision run, converted to binary64 for summaries; high-precision final states retained in JSON.'))
    dump(target/'arrays.json',dict(files=arrays,total_bytes=sum(r['bytes'] for r in arrays)))
    dump(target/'original-summaries.json',dict(files=originals))
    dump(target/'environment.json',dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__,mpmath=mpmath.__version__,matplotlib=matplotlib.__version__,threads=4,fastmath=False))
    dump(target/'workflow.json',dict(completed=['576-case pilot','40 cases continued for 500,000 steps','same 40 cases continued for another 5,000,000 steps','four full spectra','four arbitrary-precision runs','24 state controls','model and actual-state checks','40 recurrence checks'],
        development_failures=[dict(stage='Initial pilot compilation',outcome='No simulations ran; builtin abs(array) was unsupported by Numba and was replaced with np.abs(array).'),
        dict(stage='spectra01',outcome='Case 43 completed, then the next input file was unavailable while the parent extension was running. The unchanged driver was rerun to completion as spectra02 after all parents finished.'),
        dict(stage='recurrence01',outcome='PCA variance normalization underflowed for a near-zero state. Rescaling singular values fixed this diagnostic; recurrence02 is retained.'),
        dict(stage='controls01',outcome='Ordinary Sinkhorn is singular and the trial tangent could become exactly zero. Floored log norms were not meaningful growth estimates. Replayed controls02 suppresses growth reporting for all four ordinary-update controls; state results are unchanged.')],
        scope='No infinite-time chaos, basin or historical-novelty certificate. Full pilot and continuation histories were not saved; retained windows are explicitly named.'))
    print(json.dumps(dict(files=len(list(target.rglob('*'))),external_arrays=len(arrays),external_bytes=sum(r['bytes'] for r in arrays))))

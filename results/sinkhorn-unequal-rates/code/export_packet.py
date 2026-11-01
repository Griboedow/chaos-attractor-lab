"""Retain exact sources and compact JSON summaries; keep all trajectory arrays external."""
from pathlib import Path
import argparse,hashlib,json,platform,shutil
import numpy,scipy,numba,mpmath,matplotlib

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=False);work=Path.cwd();datasets=[]
    code=['model.py','model_dual.py','base_survey.py','base_followup.py','survey_rectangular.py','follow_rectangular.py','extend_rectangular.py','survey_dual.py','follow_dual.py','survey_balanced.py','follow_balanced.py','extend_selected.py','spectrum_dual.py','precision_dual.py','controls_balanced.py','check_schur.py','check_actual_dual.py','analyze_runs.py','preview_balanced.py','preview_dual.py','prepare_presentation.py','export_packet.py']
    (out/'code').mkdir()
    for name in code:shutil.copyfile(work/name,out/'code'/name)
    stages=['pilot01','followup01','rectangular-long01','dual01','dual-follow01','dual-long01','balanced01','balanced-follow01','balanced-long01','selected-100m','controls-balanced01']
    for stage in stages:
        for source in sorted((work/stage).glob('*.json')):
            data=json.loads(source.read_text());runs=data.pop('runs');stem=source.stem;parts=[]
            for k in range(0,len(runs),64 if stage in ('pilot01','dual01','balanced01') else 8):
                batch=64 if stage in ('pilot01','dual01','balanced01') else 8;rel=f'{stage}/{stem}-rows{k:04d}.json';write(out/rel,dict(runs=runs[k:k+batch]));parts.append(rel)
            index=dict(metadata=data,parts=parts,runs=len(runs),original_path=str(source),original_bytes=source.stat().st_size,original_sha256=digest(source));rel=f'{stage}/{stem}-index.json';write(out/rel,index);datasets.append(dict(original=f'{stage}/{source.name}',index=rel,original_sha256=digest(source)))
    for name in ['actual-checks.json','schur-check.json','spectra-balanced-1m.json','spectra-rectangular-6m.json','analysis.json','preview-balanced-1m.json','preview-dual-1m.json','spectra-selected-106m.json']:
        shutil.copyfile(work/name,out/name)
    for folder in ('precision-balanced01','precision-rectangular01','precision-balanced-6m'):
        (out/folder).mkdir()
        for source in sorted((work/folder).glob('*.json')):shutil.copyfile(source,out/folder/source.name)
    write(out/'datasets.json',dict(datasets=datasets,note='JSON values are preserved; long source documents are split into bounded-size parts. Original byte hashes preserve continuation lineage.'))
    arrays=[]
    for folder in stages+['precision-balanced01','precision-rectangular01','precision-balanced-6m']:
        for source in sorted((work/folder).glob('*')):
            if source.suffix not in ('.npz','.npy'):continue
            a=numpy.load(source)
            content=[dict(key=k,shape=list(a[k].shape),dtype=str(a[k].dtype)) for k in a.files] if source.suffix=='.npz' else [dict(shape=list(a.shape),dtype=str(a.dtype))]
            arrays.append(dict(path=str(source),bytes=source.stat().st_size,sha256=digest(source),arrays=content))
    write(out/'arrays.json',dict(files=arrays,total_bytes=sum(r['bytes'] for r in arrays),scope='External numerical arrays. Each tail filename names its actual retained final window. Earlier full histories were not stored; no hidden thinning. NaN rows denote windows not reached after escape or numerical capture.'))
    write(out/'environment.json',dict(python=platform.python_version(),platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__,mpmath=mpmath.__version__,matplotlib=matplotlib.__version__,numba_threads=4,selected_long_run_threads=2,precision_workers=2,blas_threads=1))
    print(json.dumps(dict(datasets=len(datasets),files=len(list(out.rglob('*.*'))),external_arrays=len(arrays),external_bytes=sum(r['bytes'] for r in arrays))))

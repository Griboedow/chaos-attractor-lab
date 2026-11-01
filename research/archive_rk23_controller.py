"""Freeze all declared controller screens, selected checks and scientific renders."""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

def arrays(path):
    with np.load(path) as f:return {k:f[k] for k in f.files}

def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);copies=[];omissions=[];splits=[]
    folders={
        'pi/scout':'adaptive-pi-scout01','pi/follow':'adaptive-pi-follow01','pi/gallery':'adaptive-pi-gallery01',
        'complex/scout':'rk23-controller-scout01','complex/follow':'rk23-controller-follow01','complex/gallery':'rk23-controller-gallery01',
        'complex/upstream-initial':'rk23-controller-upstream01','complex/upstream-control':'rk23-controller-upstream-control02','complex/upstream-analysis':'rk23-controller-upstream-analysis01',
        'oscillator/scout':'rk23-oscillator-scout01','oscillator/follow':'rk23-oscillator-follow01','oscillator/gallery':'rk23-oscillator-gallery01',
        'oscillator/five-million':'rk23-oscillator-long01','oscillator/precision':'rk23-oscillator-precision01','oscillator/cohort':'rk23-oscillator-cohort01',
        'oscillator/upstream-failed-analysis':'rk23-oscillator-upstream01','oscillator/upstream-failed-input':'rk23-oscillator-upstream02','oscillator/upstream':'rk23-oscillator-upstream03',
        'oscillator/approach':'rk23-oscillator-approach01','oscillator/render':'rk23-oscillator-render01','oscillator/evidence-initial':'rk23-oscillator-evidence01','oscillator/evidence':'rk23-oscillator-evidence02',
        'oscillator/hyper-five-million':'rk23-oscillator-hyper-long01','oscillator/hyper-precision':'rk23-oscillator-hyper-precision01',
        'oscillator/hyper-upstream':'rk23-oscillator-hyper-upstream01','oscillator/hyper-gallery':'rk23-oscillator-hyper-gallery01',
        'oscillator/error-norm-comparison':'rk23-error-norm-comparison01',
    }
    for dest,name in folders.items():
        source=args.external/name;targetdir=out/dest;targetdir.mkdir(parents=True)
        # Failed runs intentionally have failure metadata rather than a report.
        assert (source/'report.json').exists() or (source/'analysis.json').exists() or (source/'failure.json').exists(),source
        joined=arrays(source/'follow.npz') if dest in ['oscillator/follow','oscillator/five-million','oscillator/hyper-five-million'] else None
        for p in sorted(source.iterdir()):
            if not p.is_file():continue
            record=dict(source=str(p),bytes=p.stat().st_size,sha256=sha(p));reason=None
            if joined is not None and p.name.startswith('block-'):
                d=arrays(p);indices=[list(joined['ids']).index(i) for i in d['ids']]
                assert all(np.array_equal(joined[k][indices],v,equal_nan=True) for k,v in d.items())
                reason='Every array equals the corresponding rows of the archived aggregate.'
            elif dest=='pi/follow' and p.name=='follow.npz':
                d=arrays(p);count=len(d['ids']);assert count==1089
                for start in range(0,count,300):
                    stop=min(start+300,count);target=targetdir/f'follow-{start:04}-{stop:04}.npz'
                    np.savez_compressed(target,**{k:v[start:stop] for k,v in d.items()})
                    back=arrays(target);assert all(np.array_equal(v[start:stop],back[k],equal_nan=True) for k,v in d.items())
                    splits.append(dict(source=str(p),source_sha256=record['sha256'],start=start,stop=stop,path=target.relative_to(out).as_posix(),sha256=sha(target)))
                reason='Losslessly split by rows to stay below the 100 MB individual-file limit; every array compared.'
            elif dest in ['pi/gallery','complex/gallery'] and p.name.startswith('case-') and p.suffix=='.npz':
                reason='Non-promoted geometry cloud retained locally. All numerical outcomes, initial states, figure cameras, generation source and rendered gallery pages are archived; this file can be regenerated.'
            elif dest=='oscillator/render' and p.suffix=='.npy':
                reason='Eight-million-point raw render input retained locally. Exact histograms, consecutive animation states, initial data, source and trajectory hash are archived.'
            elif dest=='oscillator/upstream-failed-analysis' and p.name=='case-381.npz':
                assert sha(args.external/'rk23-oscillator-upstream03'/p.name)==record['sha256']
                reason='Byte-identical trajectory copied into the completed upstream/ run; no actual solver recomputation.'
            if reason:omissions.append(dict(**record,reason=reason));continue
            target=targetdir/p.name;shutil.copy2(p,target);assert sha(target)==record['sha256'];copies.append(dict(**record,path=target.relative_to(out).as_posix()))
    source_dir=out/'source';source_dir.mkdir()
    for path in [args.external/'adaptive-source/sunadaptcontroller_soderlind.c',args.external/'rk23-structure-audit.json']:
        target=source_dir/path.name;shutil.copy2(path,target);copies.append(dict(source=str(path),path=target.relative_to(out).as_posix(),sha256=sha(path),bytes=path.stat().st_size))
    sources=[]
    descriptions=[
        ('aves-griffiths-higham-1997.pdf','https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P28.pdf','Abstract, introduction, fixed-point/error-control discussion and conclusions; not a full-page-by-page review.'),
        ('hall-higham-1988.pdf','https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P1.pdf','Coupled polynomial recurrence, transformed norm, stability result and Section 3 practical norm comparison.'),
        ('higham-humphries-wain-2000.pdf','https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P35.pdf','Abstract and stable-equilibrium/EPS discussion in Section 2, plus reference chain; not an exhaustive reading.'),
        ('higham-hall-1990.pdf','https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P6.pdf','Visually inspected the first two scanned pages. Text extraction was empty. No claim of full-text review.'),
    ]
    for name,url,scope in descriptions:
        path=args.external/'adaptive-source'/name;sources.append(dict(title_file=name,url=url,sha256=sha(path),bytes=path.stat().st_size,read_scope=scope,redistributed=False))
    (source_dir/'primary-inspection.json').write_text(json.dumps(dict(date='2026-11-01',sources=sources,additional_web_sources='See docs/rk23-priority.md for primary web sources, bounded reading scopes and access gaps.',sundials_git_blob_sha='11ab153a9b8762c988eeff3d8041ac8292e07460',sundials_scope='Controller-source comparison only; no complete SUNDIALS solver replay.'),indent=2),encoding='utf-8')
    assets=[]
    for folder,name in [('render','rk23-oscillator-geometry.png'),('render','rk23-oscillator-rotation.gif'),('evidence','rk23-oscillator-evidence.png'),('evidence','rk23-oscillator-approach.png')]:
        source=out/'oscillator'/folder/name;target=args.assets/name;assert not target.exists();shutil.copy2(source,target)
        assets.append(dict(path=name,archive_path=source.relative_to(out).as_posix(),sha256=sha(target),bytes=target.stat().st_size))
    (out/'assets.json').write_text(json.dumps(dict(files=assets),indent=2),encoding='utf-8')
    (out/'copy-provenance.json').write_text(json.dumps(dict(copies=copies,omitted_local_files=omissions,lossless_splits=splits),indent=2),encoding='utf-8')
    (out/'README.txt').write_text('Adaptive numerical-controller study, 2026-11-01.\nAll 8198 Heun/Euler PI settings, 1027 complex RK23 settings and 1027 real oscillator settings are retained. All scout flags received the declared follow-up, including periodic and failed outcomes.\nPI is a declared custom guard, not a complete SUNDIALS solver. Complex RK23 has exact phase symmetry. The original-coordinate real oscillator is the main visual study. The physical ODE is stable.\nArchive includes six five-million-step oscillator runs, twelve independent high-precision runs, actual SciPy trajectories, 189 start/parameter controls, 65536 volume trajectories, exact portrait histograms, animation data and frozen sources. See docs/rk23-oscillator.md and docs/rk23-priority.md for claims and limits.\nLarge raw portraits and non-promoted PI/complex gallery clouds stay local with hashes and reproduction information. No outcome is removed because it is regular or inconvenient. Failed drivers and initial figures are retained.\n',encoding='utf-8')
    for name in ['archive_rk23_controller.py','verify_rk23_controller.py','audit_rk23_structure.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    files=[dict(path=p.relative_to(out).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.rglob('*')) if p.is_file()]
    assert max(f['bytes'] for f in files)<100000000
    manifest=dict(count=len(files),bytes=sum(f['bytes'] for f in files),files=files);(out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps({k:manifest[k] for k in ['count','bytes']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--external',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--assets',type=Path,required=True);main(p.parse_args())

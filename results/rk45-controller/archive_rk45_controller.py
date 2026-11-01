"""Freeze every declared RK45 screen outcome and the selected checks."""
import argparse,json,shutil
from pathlib import Path
import numpy as np
from verify_rk23_controller import read,sha


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);copies=[];omissions=[]
    folders={
        'oscillator/scout':'rk45-oscillator-scout01','oscillator/follow':'rk45-oscillator-follow01',
        'oscillator/gallery':'rk45-oscillator-gallery01','oscillator/angles':'rk45-oscillator-angles01','oscillator/hyper-gallery':'rk45-hyper-gallery01',
        'oscillator/five-million':'rk45-oscillator-long01','oscillator/precision':'rk45-oscillator-precision01',
        'oscillator/upstream':'rk45-oscillator-upstream01','oscillator/cohort':'rk45-oscillator-cohort01',
        'oscillator/approach':'rk45-oscillator-approach01','oscillator/render':'rk45-oscillator-render01',
        'oscillator/evidence-initial':'rk45-oscillator-evidence01','oscillator/evidence':'rk45-oscillator-evidence02','oscillator/branches':'rk45-oscillator-branches01','oscillator/phase':'rk45-oscillator-phase01',
        'oscillator/error-norm-comparison':'rk45-oscillator-norms01','oscillator/damping':'rk45-rosette-damping01','oscillator/library-clouds':'rk45-library-clouds01',
        'source/structure':'rk45-symbolic','source/priority':'rk45-priority','source/driver-drafts':'rk45-driver-drafts',
    }
    for dest,name in folders.items():
        source=args.external/name;targetdir=out/dest;targetdir.mkdir(parents=True)
        assert any((source/n).exists() for n in ['report.json','analysis.json','failure.json','higham-hall-1990-full-review.json']),source
        joined=read(source/'follow.npz') if dest in ['oscillator/follow','oscillator/five-million'] else read(source/'screen.npz') if dest=='oscillator/scout' else read(source/'scan.npz') if dest=='oscillator/damping' else None
        for p in sorted(source.iterdir()):
            if not p.is_file():continue
            record=dict(source=str(p),bytes=p.stat().st_size,sha256=sha(p));reason=None
            if joined is not None and p.name.startswith('block-'):
                d=read(p)
                if 'ids' in d:indices=[list(joined['ids']).index(i) for i in d['ids']]
                else:
                    offset=int(p.stem.split('-')[1]);indices=list(range(offset,offset+len(d['params'])))
                assert all(np.array_equal(joined[k][indices],v,equal_nan=v.dtype.kind in 'fc') for k,v in d.items())
                reason='Every array equals the corresponding rows of the archived aggregate.'
            elif dest=='oscillator/render' and p.suffix=='.npy':
                reason='Eight-million-state raw render input retained locally. Exact histograms, consecutive animation states, initial data, source and trajectory hash are archived.'
            if reason:omissions.append(dict(**record,reason=reason));continue
            target=targetdir/p.name;shutil.copy2(p,target);assert sha(target)==record['sha256'];copies.append(dict(**record,path=target.relative_to(out).as_posix()))
    assets=[]
    for folder,name in [('render','rk45-oscillator-geometry.png'),('render','rk45-oscillator-rotation.gif'),('evidence','rk45-oscillator-evidence.png'),('evidence','rk45-oscillator-approach.png'),('library-clouds','rk45-library-clouds.png')]:
        source=out/'oscillator'/folder/name;target=args.assets/name;assert not target.exists();shutil.copy2(source,target)
        assets.append(dict(path=name,archive_path=source.relative_to(out).as_posix(),sha256=sha(target),bytes=target.stat().st_size))
    (out/'assets.json').write_text(json.dumps(dict(files=assets),indent=2),encoding='utf-8')
    (out/'copy-provenance.json').write_text(json.dumps(dict(copies=copies,omitted_local_files=omissions),indent=2),encoding='utf-8')
    (out/'README.txt').write_text('Standard Dormand-Prince RK45 controller study, 2026-11-01.\nAll 2051 scout settings and all 1061 flagged follow-ups are retained, including failed, regular and weak outcomes. Six selected five-million-step checks; eight independent 70/90-digit checks; four unmodified SciPy trajectories; 252 start/parameter controls; 65536 volume trajectories; 514 declared damping-slice settings.\nThe physical linear oscillator is stable. D335 and D448 are numerical-controller candidates; names describe appearance only. Exact solution-coordinate bounds do not certify a full-state attractor, asymptotic chaos or historical uniqueness.\nAll 24 gallery clouds, exact portrait histograms, consecutive animation states, rigid cameras, source snapshots and algebra certificates are public. Two raw eight-million-state portraits remain local with hashes and regeneration information. Duplicate scout/follow blocks were compared with their public aggregates before omission. No prior-paper PDFs are redistributed.\nSee docs/rk45-oscillator.md, docs/rk45-controller-structure.md and docs/rk45-priority.md.\n',encoding='utf-8')
    for n in ['archive_rk45_controller.py','verify_rk45_controller.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    files=[dict(path=p.relative_to(out).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.rglob('*')) if p.is_file()]
    assert max(f['bytes'] for f in files)<100000000
    manifest=dict(count=len(files),bytes=sum(f['bytes'] for f in files),files=files)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps({k:manifest[k] for k in ['count','bytes']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--external',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--assets',type=Path,required=True);main(p.parse_args())

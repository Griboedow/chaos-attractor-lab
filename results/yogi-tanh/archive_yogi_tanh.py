"""Freeze smooth-Yogi evidence; split the large gallery without losing a single state."""
import argparse,json,shutil
from pathlib import Path
import numpy as np
from verify_optimizer_followups import sha,arrays


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);copies=[];omitted=[];splits=[]
    folders={
        'scout':'yogi-tanh','follow':'yogi-tanh-follow','five-million':'yogi-tanh-five-million',
        'precision':'yogi-tanh-precision','cohorts':'yogi-tanh-cohorts','gallery':'yogi-tanh-gallery',
        'preview':'yogi-tanh-preview','render-initial':'yogi-tanh-render','render':'yogi-tanh-render02',
        'evidence-initial':'yogi-tanh-evidence','evidence':'yogi-tanh-evidence02','source':'yogi-tanh-source',
        'upstream-failed-graph':'yogi-tanh-upstream','upstream-failed-hyper':'yogi-tanh-upstream02',
        'upstream':'yogi-tanh-upstream03','upstream-strong':'yogi-tanh-upstream-strong',
        'upstream-analysis-failed':'yogi-tanh-upstream-analysis','upstream-analysis':'yogi-tanh-upstream-analysis02',
        'upstream-strong-analysis':'yogi-tanh-upstream-strong-analysis',
        'structure-failed-1':'yogi-tanh-structure','structure-failed-2':'yogi-tanh-structure02',
        'structure-failed-3':'yogi-tanh-structure03','structure':'yogi-tanh-structure04',
    }
    for dest,source_name in folders.items():
        source=args.external/source_name;targetdir=out/dest;targetdir.mkdir()
        joined=arrays(source/('scout.npz' if dest=='scout' else 'histories.npz')) if dest in ['scout','follow'] else None
        for p in sorted(source.iterdir()):
            if not p.is_file():continue
            record=dict(source=str(p),bytes=p.stat().st_size,sha256=sha(p))
            reason=None
            if joined is not None and p.name.startswith('block-'):
                block=arrays(p);indices=[list(joined['ids']).index(i) for i in block['ids']]
                assert all(np.array_equal(joined[k][indices],v,equal_nan=True) for k,v in block.items())
                reason='Exact duplicate rows in the archived aggregate; every array compared before omission.'
            elif p.suffix=='.npy' and dest in ['preview','render','render-initial']:
                reason='Large render input retained locally with this hash; regeneration source, initial state and protocol are archived. Publication has exact density counts and consecutive animation states.'
                if dest=='render':assert sha(args.external/'yogi-tanh-render'/p.name)==record['sha256']
            elif dest=='render-initial' and p.name.endswith('-animation.npz'):
                other=args.external/'yogi-tanh-render02'/p.name
                assert all(np.array_equal(v,arrays(other)[k]) for k,v in arrays(p).items())
                reason='Exact duplicate animation arrays archived in render/; compared elementwise.'
            elif dest=='gallery' and p.name=='histories.npz':
                d=arrays(p);assert d['states'].shape==(36,200000,3)
                for i,ident in enumerate(d['ids']):
                    target=targetdir/f'case-{int(ident)}.npz';np.savez_compressed(target,**{k:v[i] for k,v in d.items()})
                    reconstructed=arrays(target);assert all(np.array_equal(v[i],reconstructed[k]) for k,v in d.items())
                    splits.append(dict(source=str(p),source_sha256=record['sha256'],index=i,id=int(ident),destination=target.relative_to(out).as_posix(),sha256=sha(target)))
                reason='Losslessly split into 36 per-case NPZ files to respect the 100 MB GitHub file limit. Every array compared elementwise.'
            if reason:omitted.append(dict(**record,reason=reason));continue
            target=targetdir/p.name;shutil.copy2(p,target);assert sha(target)==record['sha256']
            copies.append(dict(source=str(p),destination=target.relative_to(out).as_posix(),sha256=record['sha256']))
    assets=[]
    for folder,name in [('render','yogi-tanh-geometry.png'),('render','yogi-tanh-rotation.gif'),('evidence','yogi-tanh-evidence.png')]:
        p=out/folder/name;target=args.assets/name;assert not target.exists();shutil.copy2(p,target)
        assets.append(dict(path=name,archive_path=folder+'/'+name,bytes=target.stat().st_size,sha256=sha(target)))
    (out/'assets.json').write_text(json.dumps(dict(files=assets),indent=2)+'\n')
    (out/'copy-provenance.json').write_text(json.dumps(dict(files=copies,omitted_local_inputs=omitted,lossless_gallery_split=splits),indent=2)+'\n')
    (out/'README.txt').write_text('Existing TensorFlow Addons tanh Yogi: scalar-quadratic geometry and finite expansion evidence.\nAll 8256 scouts, all 488 selected followups, all 420 cohorts, six precision runs and ten actual library runs are retained, including regular and failed outcomes.\nGallery clouds are losslessly split per case. Large raw preview/render arrays remain local with hashes; sources, protocols, portrait histograms and consecutive animation arrays are retained. Duplicate block arrays were checked against archived aggregates.\nInitial rendering and failed control drivers are retained separately. The gallery screen reflects its stored camera second coordinate because density already flips the vertical image and the gallery uses origin=lower; this is a screen orientation, not a state transformation or nonlinear deformation.\nSee docs/yogi-tanh.md and docs/yogi-tanh-structure.md. No new-optimizer or historical-uniqueness claim.\n')
    for name in ['archive_yogi_tanh.py','verify_yogi_tanh.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    files=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    assert max(r['bytes'] for r in files)<100_000_000
    (out/'manifest.json').write_text(json.dumps(dict(files=files,count=len(files),bytes=sum(r['bytes'] for r in files)),indent=2)+'\n')
    print('Frozen',len(files),'files;',sum(r['bytes'] for r in files),'bytes',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['external','output','assets']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())

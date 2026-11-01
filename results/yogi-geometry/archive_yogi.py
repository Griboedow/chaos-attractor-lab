"""Freeze completed Yogi research without discarding unsuccessful experiments."""
import argparse,hashlib,json,shutil
from pathlib import Path


FOLDERS=['scout','follow-million','gallery','spectra-million','spectra-five-million','cold-five-million',
    'precision','precision-cold','precision-square','cohorts','previews','render-input','render',
    'boundary-audit','structure','upstream','upstream02','upstream-analysis','upstream-analysis02','evidence','evidence02','source']


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);external=[];copied=[]
    for name in FOLDERS:
        source=args.external/('yogi-'+name);assert source.is_dir(),source
        for p in sorted(source.iterdir()):
            if not p.is_file():continue
            relative=Path(name)/p.name
            if (name=='source' and p.suffix in ['.pdf','.txt']) or p.suffix=='.npy':
                external.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),
                    reason='Published full text retained locally; cite original source.' if name=='source' else 'Large full portrait trajectory; published histograms and animation clouds retained.'))
                continue
            target=out/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
            assert sha(p)==sha(target);copied.append(dict(source=str(p),destination=relative.as_posix(),sha256=sha(p)))
    (out/'external-inventory.json').write_text(json.dumps(dict(files=external),indent=2)+'\n')
    (out/'copy-provenance.json').write_text(json.dumps(dict(files=copied),indent=2)+'\n')
    reviews=dict(evidence='The first evidence composition has insufficient footer margin; evidence02 fixes layout without changing data.',
        upstream='The first reference assumed zero first moment; installedOptax0.2.8 initializes both memories to1e-6. The initial assertion failed before any long experiment. upstream02 is accepted.',
        upstream_analysis='The original fixed-line tangent frame was coordinate-aligned. upstream-analysis02 uses a generic frame and exact linear controls; preserve both.',
        initial_switch='precision records exact promotion of stored coordinates. precision-square explicitly defines v0=x0*x0 at working precision. They are different initial specifications; the published sunburst anchor uses zero-memory7590.',
        boundary_audit='First-step branch comparison was an inline calculation. Its exact rational switching-margin conclusion is independently reproducible by structure/analyze_yogi_structure.py.',
        render_input='Endpoints were assembled by exact copying. provenance.json and the verifier identify and check every copied source row.')
    (out/'review.json').write_text(json.dumps(reviews,indent=2)+'\n')
    (out/'README.txt').write_text('Frozen Yogi geometry study, 1 November2026.\n'
        'Known Zaheer2018Algorithm2 on theta^2/2; numerical candidates, no historical uniqueness or asymptotic chaos theorem.\n'
        'All8202scouts,268follow-ups,24gallery clouds,selected spectra,472cohorts,12precision runs and9actual-library histories are retained.\n'
        'Five-million spectra: spectra-five-million contains44,7591,4214,5174; cold-five-million contains the final zero-memory sunburst7590.\n'
        'precision7591 stored-coordinate outcomes settle; precision-square7591 relation-defined outcomes remain expanding. Both are intentional records.\n'
        'upstream and upstream-analysis are failed/superseded attempts. Accepted directories are upstream02 and upstream-analysis02.\n'
        'render contains full1m-state animation clouds and8m-state portrait histograms. Full8m/2m arrays are local-only, inventoried with hashes.\n'
        'source contains metadata and inspected implementation snapshots, not redistributed published PDFs.\n'
        'Paths in executed protocols identify original local inputs. copy-provenance.json maps archived files byte-for-byte.\n'
        'review.json documents corrections. All executed scripts remain unchanged in their individual experiment folders.\n'
        'Run python research/verify_yogi_packet.py --root results/yogi-geometry from the repository root.\n')
    for name in ['archive_yogi.py','verify_yogi_packet.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    assets=args.assets;assets.mkdir(parents=True,exist_ok=True);links=[]
    for folder,name in [('render','yogi-geometry.png'),('render','yogi-rotation.gif'),('evidence02','yogi-two-directions.png')]:
        target=assets/name;assert not target.exists();shutil.copy2(out/folder/name,target)
        links.append(dict(path=name,archive_path=f'{folder}/{name}',bytes=target.stat().st_size,sha256=sha(target)))
    (out/'assets.json').write_text(json.dumps(dict(files=links),indent=2)+'\n')
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    (out/'manifest.json').write_text(json.dumps(dict(files=entries,count=len(entries),bytes=sum(v['bytes'] for v in entries)),indent=2)+'\n')
    print(json.dumps(dict(files=len(entries),bytes=sum(v['bytes'] for v in entries),external_files=len(external))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--assets',type=Path,required=True);main(p.parse_args())

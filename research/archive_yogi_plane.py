"""Freeze the Yogi grid extension and scalar boundedness controls."""
import argparse,hashlib,json,shutil
from pathlib import Path


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);copied=[];external=[]
    folders=['plane','plane-gallery','plane-five-million','plane-preview','plane-precision','plane-upstream',
             'plane-upstream-analysis','plane-arc','plane-arc02','plane-evidence']
    pairs=[(name,args.external/('yogi-'+name)) for name in folders]
    pairs += [('scalar-bounds',args.external/'optimizer-scalar-bounds'),('bound-sources',args.external/'optimizer-bound-source')]
    for folder,source in pairs:
        assert source.is_dir()
        for p in sorted(source.iterdir()):
            if not p.is_file():continue
            if p.suffix=='.npy' or (folder=='bound-sources' and p.suffix in ['.pdf','.txt']):
                external.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),
                    reason='Published full text retained locally and cited at its original source.' if folder=='bound-sources' else 'Large full portrait cloud; unthinned animation tail and complete image histograms retained.'))
                continue
            relative=Path(folder)/p.name;target=out/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
            assert sha(p)==sha(target);copied.append(dict(source=str(p),destination=relative.as_posix(),sha256=sha(p)))
    (out/'external-inventory.json').write_text(json.dumps(dict(files=external),indent=2)+'\n')
    (out/'copy-provenance.json').write_text(json.dumps(dict(files=copied),indent=2)+'\n')
    reviews=dict(parameter_plot='plane-gallery/yogi-parameter-plane.png uses an imshow extent that shifts outer cell centers by half a bin. plane-evidence replaces only this plot with correctly centered pcolormesh cells; all original data and gallery pages retained.',
        arc_render='plane-arc completed and its explicit point-count checks passed, but the symbolic vertical-envelope formula permuted camera axes. plane-arc02 corrects that formula and the image y-axis convention, reusing exactly the same full cloud. Only plane-arc02 is published.',
        proof='scalar-bounds/proof.md is the text at execution time. The current repository note clarifies its introductory wording; the theorem, constants and proof are unchanged.',
        identifiers='G-prefixed grid IDs are source-row indices within this 2048-point grid, distinct from original scout case IDs.',
        negative_outcomes='G0063 has a nonpositive leading window in the five-million autonomous run; its actual-library history differs. G1506 settles to the fixed line in actual Optax although its autonomous run remains expanding. Neither is discarded.',
        sources='Only metadata from bound-sources is redistributed. Some downloaded sources were not read in full; inspection scope is recorded individually. No historical uniqueness conclusion follows.')
    (out/'review.json').write_text(json.dumps(reviews,indent=2)+'\n')
    (out/'README.txt').write_text('Yogi parameter-grid extension and scalar boundedness controls, 1 November 2026.\n'
        'All 2048 independent cold starts, 24 gallery clouds, six five-million runs, four 70/90-digit checks and six actual Optax histories retained.\n'
        'The selected feathered arc is G0670. G0125 is the stronger two-direction numerical candidate. Grid identifiers are not original scout identifiers.\n'
        'plane-arc02 and plane-evidence contain accepted figures. Earlier compositions and all negative outcomes remain archived; see review.json.\n'
        'The scalar boundedness proof applies to the stated scalar recurrences; exact algebra and finite stress controls do not replace the proof.\n'
        'All executed sources are frozen byte-for-byte. Source paths identify original local inputs; copy-provenance.json maps archived copies.\n'
        'Full publication PDFs and large .npy portrait clouds remain local, with checksums in external-inventory.json.\n'
        'Run python research/verify_yogi_plane.py --root results/yogi-plane from the repository root.\n')
    links=[]
    for folder,name in [('plane-arc02','yogi-feathered-arc.png'),('plane-arc02','yogi-feathered-arc.gif'),
                        ('plane-evidence','yogi-parameter-plane.png'),('plane-evidence','yogi-grid-directions.png')]:
        target=args.assets/name;assert not target.exists();shutil.copy2(out/folder/name,target)
        links.append(dict(path=name,archive_path=f'{folder}/{name}',bytes=target.stat().st_size,sha256=sha(target)))
    (out/'assets.json').write_text(json.dumps(dict(files=links),indent=2)+'\n')
    for name in ['archive_yogi_plane.py','verify_yogi_plane.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    (out/'manifest.json').write_text(json.dumps(dict(files=rows,count=len(rows),bytes=sum(r['bytes'] for r in rows)),indent=2)+'\n')
    print('Archived',len(rows),'files',sum(r['bytes'] for r in rows),'bytes',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['external','output','assets']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())

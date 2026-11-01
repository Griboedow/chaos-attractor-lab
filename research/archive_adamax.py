"""Freeze Adamax evidence, explicit failures and actual-state visual assets."""
import argparse,hashlib,importlib.metadata,json,platform,shutil,sys
from pathlib import Path


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    root=Path(__file__).resolve().parents[1];out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    mapping=[('scout','adamax-scout'),('gallery-rejected','adamax-gallery'),('gallery','adamax-gallery02'),
        ('long-million','adamax-long'),('long-four-million','adamax-four-million'),('precision-launch-failed','adamax-precision'),
        ('precision','adamax-precision02'),('cohorts','adamax-cohorts'),('upstream-short','adamax-upstream'),
        ('upstream-long','adamax-upstream-long'),('upstream-analysis-failed','adamax-upstream-analysis'),
        ('upstream-analysis','adamax-upstream-analysis02'),('structure','adamax-structure'),('preview','adamax-previews'),
        ('render','adamax-render'),('approach','adamax-approach'),('approach-layout-revision','adamax-approach-render02'),
        ('approach-low-contrast','adamax-approach-render03'),('approach-final','adamax-approach-render04'),('sheets','adamax-sheets'),
        ('outliers','adamax-outliers'),('outlier-analysis-failed','adamax-outliers-analysis02'),
        ('outlier-analysis','adamax-outliers-analysis03'),('cycle-symmetry','adamax-cycle-symmetry')]
    for name in ['adamax-four-million','adamax-precision02','adamax-cohorts','adamax-upstream-analysis02',
                 'adamax-approach','adamax-outliers-analysis03','adamax-cycle-symmetry']:
        assert (args.external/name/'analysis.json').is_file(),name
    inventory=[]
    for label,name in mapping:
        target=out/label;target.mkdir();source=args.external/name
        for path in sorted(source.iterdir()):
            if not path.is_file():continue
            item=dict(run=label,local_name=path.name,bytes=path.stat().st_size,sha256=digest(path))
            excluded=path.suffix=='.npy' or (label=='gallery-rejected' and path.suffix=='.npz')
            if excluded:
                item.update(archived=None,mode='Local only: full render/preview arrays or scientifically rejected gallery arrays. Exact generation metadata and hashes retained.')
            else:
                dest=target/path.name;shutil.copyfile(path,dest);item.update(archived=dest.relative_to(out).as_posix(),mode='byte-for-byte')
            inventory.append(item)
    sources=[];target=out/'sources';target.mkdir()
    for path in sorted((args.external/'adamax-source').iterdir()):
        if not path.is_file():continue
        sources.append(dict(name=path.name,bytes=path.stat().st_size,sha256=digest(path),included=path.suffix=='.json'))
        if path.suffix=='.json':shutil.copyfile(path,target/path.name)
    write(out/'source-inventory.json',dict(files=sources,scope='Primary PDFs, extracted text and inspection images remain local. Only retrieval metadata is redistributed; citations and inspected sections are in docs/adamax-priority.md.'))
    write(out/'local-inventory.json',inventory)
    assets=[]
    for folder,name in [('adamax-render','adamax-geometry.png'),('adamax-render','adamax-rotation.gif'),
                        ('adamax-approach-render04','adamax-approach.png'),('adamax-sheets','adamax-sheets.png')]:
        dest=root/'docs/assets'/name;assert not dest.exists();shutil.copyfile(args.external/folder/name,dest)
        assets.append(dict(path=dest.relative_to(root).as_posix(),bytes=dest.stat().st_size,sha256=digest(dest)))
    write(out/'asset-manifest.json',dict(paths_relative_to='repository root',files=assets))
    write(out/'environment.json',dict(python=sys.version,platform=platform.platform(),packages={name:importlib.metadata.version(name) for name in ['numpy','scipy','numba','matplotlib','Pillow']},
        upstream_environment='Actual PyTorch version, git version and source hash are in upstream-long/protocol.json. Arbitrary-precision checks used a separate environment with mpmath1.3.0.'))
    write(out/'review.json',dict(scope='Finite numerical candidates and exact elementary structural arguments; no certified chaotic attractor or historical uniqueness.',
        visual_QA=['2800x1680 dense portrait','72-frame1400x840 rotation, quarter-turn frames inspected and static labels checked identical',
                   '2760x1270 final approach figure, regular outcomes marked orange','2100x980 reset-sheet decomposition'],
        failures=['Gallery01completed but rejected forNaNfeatureduplicates; review.jsonsupersedesincorrectincomplete.json.',
                  'Precision01wrongfilename failed before calculation.',
                  'Upstreamanalysis01Numbaabs-array compilation failed before calculation.',
                  'Outlier integrations completed; duplicate nearest-neighbor postprocessing was stopped and replaced by exact global deduplication.',
                  'Outlieranalysis02asserted bitwise repeat136 as least period; failed. Corrected analysis03closes least-period68 under declared tolerance.',
                  'Approachimage01crowdedrowspacing; revision02fixedlayout; revision03regular-outcome colors too dim; revision04improves contrast and uses documented point markers.'],
        negative_results=['Case1103becomesperiod4.', 'Case239hasattractingperiod4coexistence.', 'Case4099has two sign-related stable68-cycles.',
                          'Scoutflags3691and3979haveexactzero(x,m)tailsdespitepositivetangentgrowth.'],
        cohort_overlap='The weak-volume4099record is Sobolstart87, already among436far endpoints.437continuationsrepresent436distinctoriginalstarts.'))
    (out/'README.txt').write_text('Frozen scalar-quadratic Adamax study.\nAll accepted scientific NPZ arrays are copied byte-for-byte, including negatives and full gallery data.\nFull8million-state portrait arrays and2million-state preview arrays remain local; their hashes and generation settings are retained.\nEach animation retains1million consecutive unthinned states. Static histograms and all frame point-count checks are included.\nRejected galleryNPZ arrays remain local and are inventoried explicitly.\nSource papers are cited, not redistributed.\nOriginal failed attempts and superseded figure versions are preserved with review notes.\nRun research/verify_adamax_packet.py from the repository root.\n',encoding='utf8')
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(r['bytes'] for r in entries),files=entries));print('Archived',len(entries),'files,',sum(r['bytes'] for r in entries),'bytes',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

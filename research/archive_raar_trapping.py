"""Freeze the strict closure run and its bounded development provenance."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
from audit_raar_trapping_batch import file_hash


def copy_tree(source,target):
    shutil.copytree(source,target,ignore=shutil.ignore_patterns('__pycache__','*.pyc','*.nbc','*.nbi'))


def main(args):
    root=Path(__file__).resolve().parents[1];src=args.base;out=args.out
    closure=json.loads((src/'raar-trapping-closure-audit01/report.json').read_text())
    assert closure['status']=='closed_accounting_verified' and closure['source_cells']==21405802 and closure['new_unprocessed_cells']==0
    out.mkdir(parents=True,exist_ok=False);(out/'candidate').mkdir()
    candidate=src/'raar-trapping-audit-candidate01'
    for name in ['protocol.json','prepare_raar_trapping_audit.py']:
        shutil.copy2(candidate/name,out/'candidate'/name)
    with (candidate/'mask.npy').open('rb') as stream, (out/'candidate/mask.npy.gz').open('wb') as destination:
        with gzip.GzipFile(filename='',fileobj=destination,mode='wb',mtime=0,compresslevel=9) as compressed:
            shutil.copyfileobj(stream,compressed)
    runs=[f'raar-trapping-audit-round2-part{i}' for i in range(4)]+[f'raar-trapping-audit-round1-part{i}' for i in range(4)]+['raar-trapping-round-comparison01','raar-trapping-closure-audit01','raar-trapping-audit-pilot01',
          'raar-interval64-validation01','raar-geometric-validation01','raar-geometric-validation02','raar-geometric-validation03',
          'raar-interval-grid-validation01','raar-interval-grid-validation02','raar-interval-grid-validation03','raar-interval-grid-validation04','raar-interval-grid-validation05','raar-interval-grid-validation06',
          'raar-trapping-render01','raar-trapping-render02']
    provenance=[]
    for name in runs:
        copy_tree(src/name,out/'runs'/name)
        for path in (src/name).rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix not in ['.pyc','.nbi','.nbc']:
                provenance.append({'path':f'runs/{name}/{path.relative_to(src/name).as_posix()}','source':str(path),'sha256':file_hash(path)})
    # Search arrays can be large and are not needed to verify the explicit set.
    # Retain all available protocol/source files from the declared proposal runs.
    proposals=['raar-tube-proposal01','raar-tube-proposal02','raar-tube-proposal03',
               'raar-trapping-grid-020-01','raar-trapping-grid-010-01','raar-trapping-grid-010-02','raar-trapping-grid-005-01',
               'raar-trapping-zono-010-01','raar-trapping-zono-010-02','raar-trapping-zono-005-01',
               'raar-trapping-bitset-control01','raar-trapping-adaptive-005-01','raar-trapping-adaptive-005-resume01','raar-trapping-adaptive-005-resume02']
    omitted=[]
    for name in proposals:
        source=src/name
        if not source.exists():
            omitted.append({'run':name,'reason':'No directory at the declared scratch location.'});continue
        (out/'proposals'/name).mkdir(parents=True)
        for path in source.iterdir():
            if not path.is_file():continue
            if path.suffix in ['.json','.jsonl','.py']:
                dest=out/'proposals'/name/path.name;shutil.copy2(path,dest)
                provenance.append({'path':dest.relative_to(out).as_posix(),'source':str(path),'sha256':file_hash(path)})
            else:
                omitted.append({'run':name,'file':path.name,'bytes':path.stat().st_size,'sha256':file_hash(path),'reason':'Proposal-only array or auxiliary file; retained locally, not required by strict replay.'})
    (out/'proposal-omissions.json').write_text(json.dumps(omitted,indent=2)+'\n',encoding='utf-8')
    shutil.copy2(src/'raar-trapping-priority.json',out/'priority.json')
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    (out/'packet.json').write_text(json.dumps({'source_cells':closure['source_cells'],'status':closure['status'],
        'candidate_storage':'Lossless deterministic gzip of mask.npy; keys.npy reconstructed from its set bits and checked against its original SHA-256.',
        'symbolic_dependency':'results/raar-symbolic-dynamics','symbolic_manifest_sha256':file_hash(root/'results/raar-symbolic-dynamics/manifest.json'),
        'scope':'Strict outward closure and exhaustive source accounting. Proposal-only arrays omitted explicitly; the candidate and all proof inputs are reconstructible. No transitivity, cloud equality or historical uniqueness assertion.'},indent=2)+'\n',encoding='utf-8')
    for name in ['region','approach']:
        shutil.copy2(src/'raar-trapping-render02'/f'raar-trapping-{name}.png',root/'docs/assets'/f'raar-trapping-{name}.png')
    print(json.dumps({'packet':str(out),'compressed_mask_bytes':(out/'candidate/mask.npy.gz').stat().st_size,'runs':len(runs),'proposals':len(proposals),'note':'Write README and priority record, then generate final manifest.'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,default=Path('results/raar-trapping'))
    main(p.parse_args())

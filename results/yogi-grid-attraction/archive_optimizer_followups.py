"""Freeze the Yogi attraction controls and the separate AdaDelta priority audit."""
import argparse,hashlib,json,shutil
from pathlib import Path


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def freeze(root,external,pairs,asset_pair,assets,description,review,proof=None):
    root.mkdir(parents=True,exist_ok=False);copied=[];omitted=[]
    for destination,source_name in pairs:
        source=external/source_name
        for p in sorted(source.iterdir()):
            if not p.is_file():continue
            if p.name=='inspected-function.py':
                omitted.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),reason='Inspected installed library implementation kept locally; its source/function hashes are recorded in analysis.json.'));continue
            relative=Path(destination)/p.name;target=root/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
            copied.append(dict(source=str(p),destination=relative.as_posix(),sha256=sha(p)));assert sha(target)==sha(p)
    if proof:shutil.copy2(proof,root/'proof-note.md')
    (root/'README.txt').write_text(description+'\n')
    (root/'review.json').write_text(json.dumps(review,indent=2)+'\n')
    (root/'copy-provenance.json').write_text(json.dumps(dict(files=copied),indent=2)+'\n')
    (root/'external-inventory.json').write_text(json.dumps(dict(files=omitted),indent=2)+'\n')
    folder,name=asset_pair;target=assets/name;assert not target.exists();shutil.copy2(root/folder/name,target)
    (root/'assets.json').write_text(json.dumps(dict(files=[dict(path=name,archive_path=f'{folder}/{name}',sha256=sha(target),bytes=target.stat().st_size)]),indent=2)+'\n')
    for name in ['archive_optimizer_followups.py','verify_optimizer_followups.py']:(root/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    entries=[dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(root.rglob('*')) if p.is_file()]
    (root/'manifest.json').write_text(json.dumps(dict(files=entries,count=len(entries),bytes=sum(r['bytes'] for r in entries)),indent=2)+'\n')
    print(root.name,len(entries),sum(r['bytes'] for r in entries),'bytes',flush=True)


def main(args):
    freeze(args.results/'yogi-grid-attraction',args.external,
        [('cohorts','yogi-grid-cohorts'),('approach','yogi-grid-approach'),('figure02','yogi-grid-approach02')],
        ('figure02','yogi-grid-approach.png'),args.assets,
        'Yogi grid attraction controls,1 November2026.236 cohort runs;32768 initial full-state points per grid candidate G0670/G0125;256 independent volume-rate checks. All data retained. Original approach figure has insufficient row-label clearance; accepted figure02 adds explicit final zooms and numerical fixed-tail markers. Run research/verify_optimizer_followups.py --results results. Finite attraction evidence, not a basin measure or chaos theorem.',
        dict(original_figure='Initial layout puts the second row heading over the first image row. Figure02 corrects spacing and adds clearly labeled final zooms with included/excluded point counts. It reuses every computed state.',
             outcomes='G0670 has1163 detected fixed tails of32768;G0125 has1968. No nonfixed endpoint is farther than.02 normalized reference extent. This is a finite sampled-support diagnostic.'))
    freeze(args.results/'adadelta-cycle-audit',args.external,
        [('replay','adadelta-cycle-audit'),('independent','adadelta-cycle-independent'),('source','adadelta-cycle-source'),('figure','adadelta-cycle-figure')],
        ('figure','adadelta-cycle-drift.png'),args.assets,
        'AdaDelta positive-epsilon cycle audit,1 November2026. Full printed recurrence checked against Bock dissertation equation5.2,Theorem5.4,Figure52. Five100-million-step numerical controls,70/90-digit5000-step replays,and actual PyTorch5000-step history retained. Exact positive-E period-two condition alpha*curvature=2 is proved in docs/adadelta-cycle-audit.md; numerical traces do not prove global convergence. Run research/verify_optimizer_followups.py --results results.',
        dict(scope='The zero-epsilon cycle family is reproduced. The positive-epsilon classification follows from all three return equations, not from a reduced fixed-memory system. No global convergence theorem or new attractor is claimed.',
             source='Source full text and page renders remain in the local optimizer-thesis-source folder; the original URL,checksum and exact inspection locations are retained.',
             underflow='The larger-E binary64 control reaches subnormal floating-point values; not an exact-real stopping result.'),Path(__file__).parents[1]/'docs/adadelta-cycle-audit.md')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['external','results','assets']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())

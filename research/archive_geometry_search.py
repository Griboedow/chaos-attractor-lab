"""Freeze compact evidence for the 2026-10-31 geometry-led screening block."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    d = args.external.resolve()
    sources = {
        'motor-pair': ROOT/'outputs/motor-geometry-20261031',
        'motor-N3': ROOT/'outputs/motor-shapes-20261031-run02',
        'straight-N3-equal': ROOT/'outputs/follower-links-20261031-equal',
        'straight-N3-stiffness': ROOT/'outputs/follower-links-20261031-stiffness',
        'straight-N3-recovery': d/'stiffness-recovery',
        'curved-N3-pilot': ROOT/'outputs/curved-follower-20261031-pilot',
        'curved-N3-broad': d/'curved-broad-run02',
        'curved-N3-check92': ROOT/'outputs/curved-follower-20261031-check92',
        'follower-N4': d/'follower-N4',
        'follower-N4-check824': d/'follower-N4-check824',
    }
    inventory = []
    for label,source in sources.items():
        target = out/label
        target.mkdir()
        for pattern in ['*.json','*.py','*.cl']:
            for path in sorted(source.glob(pattern)):
                shutil.copyfile(path,target/path.name)
        for name in ['declared-inputs.npz','complete.npz','trace.npz']:
            path = source/name
            if path.exists():
                shutil.copyfile(path,target/name)
        # All large local trajectories are indexed even where not included in Git.
        for path in sorted(source.glob('*.npz')):
            inventory.append(dict(run=label,name=path.name,bytes=path.stat().st_size,
                                  sha256=digest(path),included=(target/path.name).exists()))
    # Reconcile all 432 stiffness jobs against original or recovered raw histories.
    records = json.loads((d/'stiffness-recovery/combined-analysis.json').read_text())['runs']
    assert len(records)==432 and sorted(r['id'] for r in records)==list(range(432))
    for record in records:
        base = d/'stiffness-recovery' if record['id'] in [247,249,252] else sources['straight-N3-stiffness']
        path = base/record['file']
        assert path.is_file() and digest(path)==record['sha256'],str(path)
        with np.load(path) as data:
            for key in data.files:
                assert np.isfinite(data[key]).all(),(path,key)
    curved = np.load(out/'curved-N3-broad/complete.npz')
    quad = np.load(out/'follower-N4/complete.npz')
    assert curved['parameters'].shape==(4096,9) and not curved['retained'].any()
    assert quad['parameters'].shape==(2048,13) and np.flatnonzero(quad['retained']).tolist()==[824]
    assert np.array_equal(np.load(out/'curved-N3-pilot/complete.npz')['parameters'],curved['parameters'][:512])
    write = lambda p,x:p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')
    write(out/'local-inventory.json',inventory)
    write(out/'audit.json',dict(stiffness_histories_hash_checked=432,curved_scouts=4096,curved_retained=0,
        four_link_scouts=2048,four_link_retained=[824],pilot_is_prefix=True,
        scope='Compact published evidence. Large motor and straight-chain histories remain locally indexed, not included here.',
        failures=[dict(run='motor-shapes-20261031',reason='Empty event array shape; corrected before separate run02.'),
                  dict(run='follower-links-20261031-stiffness',reason='Disk full; 429 valid histories and three exact-job recoveries, IDs 247/249/252.'),
                  dict(run='curved-follower-20261031-broad',reason='Disk full at first observed checkpoint; separate complete run02 on D drive. No failed interval counted as negative.')]))
    entries = [dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p))
               for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(files=entries,count=len(entries),bytes=sum(p['bytes'] for p in entries)))
    print(json.dumps(dict(files=len(entries),bytes=sum(p['bytes'] for p in entries),checks='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--external',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())

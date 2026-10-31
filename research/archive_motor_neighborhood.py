"""Freeze the completed local motor survey without changing the earlier packet."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf8')


def main(args):
    out=args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    inventory=[]
    for label, name in [('scout','motor-ribbon-neighborhood'),
                        ('checks','motor-ribbon-neighborhood-checks'),
                        ('gallery','motor-ribbon-neighborhood-gallery')]:
        source=args.external/name; target=out/label; target.mkdir()
        for path in sorted(source.iterdir()):
            if path.suffix in ['.json','.py','.cl']:
                shutil.copyfile(path,target/path.name)
            elif path.suffix=='.npz':
                keep=label=='scout'
                if keep: shutil.copyfile(path,target/path.name)
                inventory.append(dict(run=label,name=path.name,bytes=path.stat().st_size,
                                      sha256=digest(path),included=keep))
    data=np.load(out/'scout/complete-or-budgeted.npz')
    assert len(data['status'])==160 and np.all(data['status']==0)
    retained=(data['lle']>.01)&(data['lle_windows'].min(axis=1)>.005)
    assert retained[:32].sum()==32 and retained[32:].sum()==84
    pairs=retained[32:].reshape(-1,2)
    assert np.all(pairs[:,0]==pairs[:,1]) and pairs.all(axis=1).sum()==42
    report=json.loads((out/'checks/analysis.json').read_text())
    assert report['completed']==report['declared']==7
    for run in report['runs']:
        path=args.external/'motor-ribbon-neighborhood-checks'/run['file']
        assert digest(path)==run['sha256']
        case=np.load(path)
        assert np.isfinite(case['states']).all()
        assert run['trace_residual_per_tau']<1e-10
        assert run['maximum_disk_radius_squared']<1+run['parameters'][2]**2
    write(out/'local-inventory.json',inventory)
    write(out/'audit.json',dict(completed=160,anchor_retained=32,anchor_declared=32,
          neighbor_lanes_retained=84,neighbor_lanes_declared=128,
          neighbor_tuples_both_retained=42,neighbor_tuples_declared=64,
          independent_spectra=7,scope='Finite-time evidence; no basin probability, open-set proof, common-topology claim, or priority claim.'))
    root=Path(__file__).resolve().parents[1]
    asset=root/'docs/assets/motor-ribbon-neighborhood.png'
    assert not asset.exists()
    shutil.copyfile(args.external/'motor-ribbon-neighborhood-gallery/gallery.png',asset)
    write(out/'asset-manifest.json',dict(paths_relative_to='repository root',files=[
        dict(path=asset.relative_to(root).as_posix(),bytes=asset.stat().st_size,sha256=digest(asset))]))
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p))
             for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries))
    print(json.dumps(dict(files=len(entries),bytes=sum(e['bytes'] for e in entries),checks='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--external',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())

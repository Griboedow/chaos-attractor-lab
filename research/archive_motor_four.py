"""Freeze the completed four-row geometry experiments and their provenance."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from PIL import Image


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    root=Path(__file__).resolve().parents[1];out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    inventory=[]
    mapping=[('scout','motor-four-asymmetry'),('checks','motor-four-checks'),('long','motor-four-long'),
        ('refined','motor-four-refined'),('source-audit','motor-four-source-audit'),('neighborhood','motor-four-neighborhood'),
        ('permutations','motor-four-permutations'),('three-row-permutations','motor-three-permutations'),
        ('gallery','motor-four-gallery'),('portrait','motor-four-portrait')]
    for label,name in mapping:
        source=args.external/name;target=out/label;target.mkdir()
        for path in sorted(source.iterdir()):
            if path.suffix in ['.json','.py','.cl']:
                shutil.copyfile(path,target/path.name)
            elif path.suffix=='.npz':
                keep=label in ['scout','neighborhood'] or (label=='long' and path.name=='case-3142.npz')
                if keep:shutil.copyfile(path,target/path.name)
                inventory.append(dict(run=label,name=path.name,bytes=path.stat().st_size,sha256=digest(path),included=keep))
    scout=np.load(out/'scout/complete-or-budgeted.npz');assert len(scout['status'])==4032 and np.all(scout['status']==0)
    neighborhood=np.load(out/'neighborhood/complete-or-budgeted.npz')
    keep=(neighborhood['status']==0)&(neighborhood['lle']>.01)&(neighborhood['lle_windows'].min(axis=1)>.005)
    assert len(keep)==160 and np.all(neighborhood['status']==0) and keep[:32].sum()==32
    assert keep[32:].sum()==82 and keep[32:].reshape(-1,2).all(axis=1).sum()==41
    for label in ['long','refined']:
        record=json.loads((out/label/'case-3142.json').read_text())
        assert digest(args.external/f'motor-four-{label}'/record['file'])==record['sha256']
        assert record['trace_residual_per_tau']<1e-10 and record['maximum_disk_radius_squared']<3.25
        assert min(w['spectrum'][0] for w in record['windows'])>.019
    perm=json.loads((out/'permutations/audit.json').read_text())
    close=[r['permutation'] for r in perm['permutations'] if r['within_three_times_identity']]
    assert close==[[0,1,2,3],[1,3,2,0],[3,0,2,1]]
    assets=[]
    for original,name in [('rotation.gif','motor-four-rotation.gif'),('poster.png','motor-four-poster.png'),('camera-audit.png','motor-four-camera-audit.png')]:
        asset=root/'docs/assets'/name;assert not asset.exists()
        shutil.copyfile(args.external/'motor-four-portrait'/original,asset);assets.append(asset)
    for path in sorted((args.external/'motor-four-gallery').glob('*.png')):
        asset=root/'docs/assets'/f'motor-four-{path.name}';assert not asset.exists()
        shutil.copyfile(path,asset);assets.append(asset)
    image=Image.open(root/'docs/assets/motor-four-rotation.gif');assert image.n_frames==72 and image.size==(800,800)
    write(out/'local-inventory.json',inventory)
    write(out/'audit.json',dict(scout_completed=4032,full_short_spectra=10,long_spectra=2,refined_spectra=1,
        neighborhood_completed=160,anchor_retained=32,neighbor_lanes_retained=82,neighbor_tuples_both_retained=41,
        empirical_close_permutations=close,animation_frames=72,
        scope='Numerical geometry and reproducibility evidence, not certified chaos, topology, attractor isotropy, or scientific priority.'))
    write(out/'asset-manifest.json',dict(paths_relative_to='repository root',files=[dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in assets]))
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries))
    print(json.dumps(dict(files=len(entries),bytes=sum(e['bytes'] for e in entries),assets=len(assets),checks='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

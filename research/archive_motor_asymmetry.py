"""Freeze completed asymmetry experiments and their visual provenance."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from PIL import Image


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);d=args.external.resolve()
    sources={'scout':'motor-asymmetry','checks':'motor-asymmetry-checks','long':'motor-asymmetry-long',
             'refined':'motor-ribbon-refined','source-audit':'motor-ribbon-audit','gallery':'motor-asymmetry-gallery02',
             'long-gallery':'motor-asymmetry-long-gallery02','camera':'motor-ribbon-camera02','animation':'motor-ribbon-animation'}
    inventory=[]
    for label,name in sources.items():
        source=d/name;target=out/label;target.mkdir()
        assert source.is_dir()
        for extension in ['json','py','cl']:
            for path in sorted(source.glob('*.'+extension)):shutil.copyfile(path,target/path.name)
        for path in sorted(source.glob('*.npz')):
            keep=label=='scout' or (label=='long' and path.name=='case-2566.npz')
            if keep:shutil.copyfile(path,target/path.name)
            inventory.append(dict(run=label,name=path.name,bytes=path.stat().st_size,sha256=digest(path),included=keep))
        if label in ['checks','long','refined']:
            analysis=json.loads((source/'analysis.json').read_text());assert analysis['declared']==analysis['completed']
            for record in analysis['runs']:
                path=source/record['file'];assert digest(path)==record['sha256']
                data=np.load(path);assert np.isfinite(data['states']).all()
                assert abs(data['log_sums'].sum()-data['integrated_trace'])/sum(np.diff([200.,data['time'][-1]]))<1e-6
    scout=np.load(out/'scout/complete-or-budgeted.npz');assert len(scout['status'])==4032 and np.all(scout['status']==0)
    primary=(scout['lle']>.005)&(scout['lle_windows'].min(axis=1)>.001)
    strict=scout['lle_windows'].min(axis=1)>.01
    assert primary.sum()==125 and strict.sum()==122
    gif=Image.open(d/'motor-ribbon-animation/rotation.gif');assert gif.n_frames==72 and gif.size==(720,720)
    frames=[]
    for i in [0,18,36,54]:gif.seek(i);frames.append(np.array(gif.convert('RGB')))
    deltas=[float(np.mean(abs(a.astype(float)-b.astype(float)))) for a,b in zip(frames,frames[1:])]
    assert min(deltas)>.2
    write=lambda p,x:p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')
    write(out/'local-inventory.json',inventory)
    write(out/'audit.json',dict(scout_lanes=4032,scout_completed=4032,preliminary=125,stricter=122,
        full_spectra_short=6,long=2,refined=1,quarter_turn_frame_mean_absolute_changes=deltas,
        included_dense_history='long/case-2566.npz',
        scope='Complete scout and selected dense history; all other local histories are hash-indexed and reproducible from saved initial states and executed sources.'))
    root=Path(__file__).resolve().parents[1]
    assets={'motor-ribbon-rotation.gif':'motor-ribbon-animation/rotation.gif',
            'motor-ribbon-poster.png':'motor-ribbon-animation/poster.png',
            'motor-ribbon-camera-audit.png':'motor-ribbon-camera02/camera-audit.png',
            'motor-asymmetry-gallery.png':'motor-asymmetry-gallery02/gallery.png'}
    asset_records=[]
    for name,source in assets.items():
        target=root/'docs/assets'/name;assert not target.exists()
        shutil.copyfile(d/source,target)
        asset_records.append(dict(path=target.relative_to(root).as_posix(),bytes=target.stat().st_size,sha256=digest(target)))
    write(out/'asset-manifest.json',dict(paths_relative_to='repository root',files=asset_records))
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries))
    print(json.dumps(dict(files=len(entries),bytes=sum(e['bytes'] for e in entries),checks='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())

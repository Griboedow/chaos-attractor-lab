"""Freeze the declared relative-smoothing mechanism control."""
import argparse,hashlib,json,shutil
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);copies=[]
    for folder,name in [('control','yogi-smoothing'),('figure','yogi-smoothing-figure'),('figure02','yogi-smoothing-figure02')]:
        source=args.external/name
        for p in sorted(source.iterdir()):
            if not p.is_file():continue
            target=out/folder/p.name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
            assert sha(p)==sha(target);copies.append(dict(source=str(p),destination=target.relative_to(out).as_posix(),sha256=sha(p)))
    (out/'copy-provenance.json').write_text(json.dumps(dict(files=copies),indent=2)+'\n')
    (out/'README.txt').write_text('Relative-width tanh smoothing of two previously checked Yogi cases.\n56declared runs, all outcomes retained. This modified equation is distinct from the existing TensorFlow Addons fixed-scale tanh option.\nfigure02 uses visible markers for detected cycles; figure retains the first density-only composition.\nNo new-optimizer, asymptotic-chaos or historical-uniqueness claim. See docs/yogi-smoothing-control.md.\n')
    asset=args.assets/'yogi-smoothing-control.png';assert not asset.exists();shutil.copy2(out/'figure02'/asset.name,asset)
    (out/'assets.json').write_text(json.dumps(dict(files=[dict(path=asset.name,archive_path='figure02/'+asset.name,sha256=sha(asset),bytes=asset.stat().st_size)]),indent=2)+'\n')
    for name in ['archive_yogi_smoothing.py','verify_yogi_smoothing.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    (out/'manifest.json').write_text(json.dumps(dict(files=rows,count=len(rows),bytes=sum(r['bytes'] for r in rows)),indent=2)+'\n');print(len(rows),'files',sum(r['bytes'] for r in rows),'bytes')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['external','output','assets']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())

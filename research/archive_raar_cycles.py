"""Freeze the complete cycle-structure follow-up without replacing the base packet."""
import argparse,json,shutil
from pathlib import Path
import numpy as np
from verify_rk23_controller import read,sha
from build_raar_cycle_figures import figures

FOLDERS={'continuation':'raar-cycles01','cycle-proof-incomplete':'raar-cycle-proof01','cycle-proof':'raar-cycle-proof02','multiplier-proof-initial':'raar-multiplier-proof01','multiplier-proof-nested':'raar-multiplier-proof02','multiplier-proof-dyadic':'raar-multiplier-proof03','multiplier-proof':'raar-multiplier-proof04','structure-audit':'raar-cycle-structure-audit01','fine-slice':'raar-cascade01','priority':'raar-cycle-priority'}


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);copies=[];omissions=[]
    for dest,name in FOLDERS.items():
        source=args.external/name;targetdir=out/dest;targetdir.mkdir();joined=read(source/'scan.npz') if dest=='fine-slice' else None
        for p in sorted(source.rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts:continue
            record=dict(source=str(p),bytes=p.stat().st_size,sha256=sha(p));reason=None
            if dest=='priority' and p.suffix!='.json':reason='Third-party publication and visual inspection rendering remain local; public metadata identifies URLs, scope and hashes.'
            elif joined is not None and p.name.startswith('block-'):
                data=read(p);indices=data['ids'];assert all(np.array_equal(joined[k][indices],v,equal_nan=v.dtype.kind in 'fc') for k,v in data.items())
                reason='Every array equals its rows in the retained scan.npz aggregate.'
            if reason:omissions.append(dict(**record,reason=reason));continue
            target=targetdir/p.relative_to(source);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);assert sha(target)==record['sha256'];copies.append(dict(**record,path=target.relative_to(out).as_posix()))
    (out/'copy-provenance.json').write_text(json.dumps(dict(copies=copies,omitted_local_files=omissions),indent=2),encoding='utf-8')
    figures(out,args.base,out/'figures')
    assets=[]
    for origin,name in [('figure-1-bifurcation.png','raar-cycle-bifurcation.png'),('figure-3-cycle-geometry.png','raar-cycle-geometry.png')]:
        source=out/'figures'/origin;target=args.assets/name;assert not target.exists();shutil.copy2(source,target);assets.append(dict(path=name,archive_path=source.relative_to(out).as_posix(),sha256=sha(target)))
    (out/'assets.json').write_text(json.dumps(dict(files=assets),indent=2),encoding='utf-8')
    (out/'README.txt').write_text('R1051 cycle-structure follow-up, 2026-11-01.\nFour numerical branches, 401 settings each; all 1604 continuation outcomes retained. A fresh-start 1025-setting fine slice retains all outcomes and its selection protocol.\nFive isolated interval certificates establish least periods 7,14,28,56 and local stability at four settings, plus an unstable least-7 cycle at the base arch. Four uniform parametric certificates establish simple minus-one multiplier events; no nondegenerate flip, infinite cascade, asymptotic chaos or historical novelty is proved.\nThe incomplete initial cycle calculation and two earlier successful multiplier calculations remain explicit. The final multiplier version adds exact dyadic endpoint serialization to the nesting check.\nThe independent interval evaluator uses the nested reflector definition. Other five geometry parameters are exact binary64 lifts; tiny interval boxes control wrapping, not physical measurement uncertainty.\nThird-party publications are linked and hashed, not redistributed. The base results/raar-circles packet remains immutable.\n',encoding='utf-8')
    for name in ['archive_raar_cycles.py','verify_raar_cycles.py','build_raar_cycle_figures.py','continue_raar_cycles.py','search_raar_circles.py','render_twist.py','verify_rk23_controller.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    files=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    manifest=dict(count=len(files),bytes=sum(r['bytes'] for r in files),base_manifest_sha256=sha(args.base/'manifest.json'),files=files)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps({k:manifest[k] for k in ['count','bytes']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--external',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--base',type=Path,default=Path('results/raar-circles'));p.add_argument('--assets',type=Path,default=Path('docs/assets'));main(p.parse_args())

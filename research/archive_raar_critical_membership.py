"""Freeze completed critical-membership evidence and explicitly retained trials."""
import argparse
import json
import shutil
from pathlib import Path
from audit_raar_coverings import sha


RUNS=['raar-critical-crossing-proposal01','raar-critical-crossing-pilot01','raar-critical-crossing-doubleton01','raar-critical-crossing-doubleton02','raar-critical-crossing-doubleton03','raar-critical-crossing-ad01','raar-critical-fold-region01','raar-critical-thin-window01','raar-critical-thin-cover01','raar-critical-thin-cover-ad01','raar-critical-thin-faces01','raar-critical-thin-faces-ad01','raar-critical-thin-faces-natural01','raar-critical-thin-fold-region01','raar-critical-thin-fold-ad01','raar-critical-membership-audit01','raar-critical-membership-audit02','raar-critical-membership-audit03','raar-critical-membership-audit04','raar-critical-membership-render01','raar-critical-membership-render02']


def main(args):
    for name in ['raar-critical-membership-audit03','raar-critical-membership-audit04']:
        assert json.loads((args.base/name/'report.json').read_text())['status']=='passed'
    args.out.mkdir(parents=True,exist_ok=False);provenance=[]
    for name in RUNS:
        folder=args.base/name;assert folder.is_dir()
        for source in sorted(folder.rglob('*')):
            if not source.is_file() or '__pycache__' in source.parts:continue
            target=args.out/'runs'/name/source.relative_to(folder);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
            provenance.append({'path':target.relative_to(args.out).as_posix(),'source':str(source),'bytes':source.stat().st_size,'sha256':sha(source)})
    (args.out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps({'runs':len(RUNS),'files':len(provenance),'bytes':sum(r['bytes'] for r in provenance),'status':'copied; add mathematical statement, priority record and manifest before publication'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

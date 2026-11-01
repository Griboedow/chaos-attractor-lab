"""Replay the declared run graph into a new external directory."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

RUNS=[
 ('pilot128','restart.py','base117',128,.006,.12,None),
 ('restart64-dt006','restart.py','base117',64,.006,96,None),
 ('restart64-dt003','restart.py','base117',64,.003,96,None),
 ('restart128-dt006','restart.py','base117',128,.006,96,None),
 ('accel128-dt006','restart_accelerated.py','base117',128,.006,96,None),
 ('field121-n64','restart_accelerated.py','base117',64,.006,96,121),
 ('field124-n64','restart_accelerated.py','base117',64,.006,96,124),
 ('rotation-seed','restart_accelerated.py','rotation',64,.004,.12,None),
 ('spectrum117-n64','tangent.py','restart64-dt006',None,None,120,24),
 ('spectrum124-n64','tangent.py','field124-n64',None,None,120,24),
 ('spectrum-rotation','tangent.py','rotation-seed',None,None,120,24),
 ('field127-from124','restart_accelerated.py','field124-n64',64,.006,120,127),
 ('field130-from124','restart_accelerated.py','field124-n64',64,.006,120,130),
 ('field124-n128-relax','restart_accelerated.py','spectrum124-n64',128,.006,24,None),
 ('spectrum124-extension','tangent_extended.py','spectrum124-n64',None,None,240,0),
 ('spectrum124-n128','tangent.py','field124-n128-relax',None,None,120,24),
 ('spectrum127-n64','tangent.py','field127-from124',None,None,120,24),
 ('field130-continuation','restart_accelerated.py','field130-from124',64,.006,120,None),
 ('field124-dt003-relax','restart_accelerated.py','spectrum124-n64',64,.003,24,None),
 ('spectrum124-dt003','tangent.py','field124-dt003-relax',None,None,120,24),
]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base117',type=Path,required=True)
    ap.add_argument('--rotation',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).resolve().parent
    seeds={'base117':args.base117.resolve(),'rotation':args.rotation.resolve()}
    for name,script,seed,n,dt,duration,extra in RUNS:
        seed_path=seeds[seed] if seed in seeds else args.output/seed
        cmd=[sys.executable,str(source/script),'--seed',str(seed_path),
             '--output',str(args.output/name),'--duration',str(duration)]
        if n is not None:
            cmd+=['--segments',str(n),'--dt',str(dt)]
            if extra is not None:cmd+=['--chi',str(extra)]
        else:
            cmd+=['--burn',str(extra)]
            if script=='tangent_extended.py':cmd+=['--resume-tangent']
        subprocess.run(cmd,cwd=args.output,check=True)
        report=json.loads((args.output/name/'report.json').read_text(encoding='utf-8'))
        if not report['success'] and name!='field130-continuation':
            raise RuntimeError(f'Unexpected incomplete integration: {name}')

if __name__=='__main__':main()

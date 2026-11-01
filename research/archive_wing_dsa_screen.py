"""Freeze finite wing-model screening and its unresolved physical scope."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);inventory=[]
    for folder in ['wing-dsa-source','wing-dsa-pilot','wing-dsa-returns','wing-dsa-long']:
        source=args.external/folder;target=out/folder;target.mkdir()
        for path in sorted(source.iterdir()):
            if path.suffix in ['.py','.cl','.json']:shutil.copyfile(path,target/path.name)
            elif path.suffix=='.npz':
                data=np.load(path);history='time' in data.files
                inventory.append(dict(folder=folder,file=path.name,bytes=path.stat().st_size,sha256=digest(path),included_in_full=not history))
                if not history:shutil.copyfile(path,target/path.name)
                else:
                    mask=data['time']>=data['time'][-1]-50
                    values={k:data[k] for k in data.files if k not in ['time','states','independent_states']}
                    values.update(time=data['time'][mask],states=data['states'][mask],first_state=data['states'][0],final_state=data['states'][-1],original_sample_count=len(data['time']))
                    if 'independent_states' in data.files:values['independent_states']=data['independent_states'][mask]
                    np.savez_compressed(target/(path.stem+'-compact.npz'),**values)
    pilot=json.loads((args.external/'wing-dsa-pilot/analysis.json').read_text())
    long=json.loads((args.external/'wing-dsa-long/analysis.json').read_text())
    returns=json.loads((args.external/'wing-dsa-returns/analysis.json').read_text())
    assert pilot['lanes']==8196 and pilot['retained_ids']==[] and len(pilot['weak_ids'])==10
    assert long['statuses']['0']==10 and long['retained_source_ids']==[] and returns['completed']==10
    write(out/'audit.json',dict(pilot_statuses=pilot['statuses'],weak_ids=pilot['weak_ids'],long_statuses=long['statuses'],long_retained=0,
        independent_return_event_counts={str(r['id']):r['events'] for r in returns['runs']},
        conclusion='No retained candidate. Independent800tau returns had too few events to classify; no periodicity claim for those cases.',
        physical_limit='Positive gain in literal2022Eq9; missing calibrated drag means no physical-angle reproduction. Does not cover the opposite-feedback-sign2023/2025 formulation.',
        extraction='All return events and final50tau at original sampling; full histories remain locally hash-indexed.'))
    write(out/'local-inventory.json',inventory);(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(r['bytes'] for r in entries),files=entries))
    print('files',len(entries),'bytes',sum(r['bytes'] for r in entries))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

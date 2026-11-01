"""Freeze the completed polymer screen, including negative and unresolved evidence."""
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
    folders=['supramolecular-source','supramolecular-pilot','supramolecular-log-recheck','supramolecular-log-long',
        'supramolecular-returns','supramolecular-slow-returns','supramolecular-budget-returns','supramolecular-2026-author']
    for folder in folders:
        source=args.external/folder;target=out/folder;target.mkdir()
        for path in sorted(source.iterdir()):
            if path.suffix in ['.json','.py','.cl'] or path.name in ['README.md','LICENCE','LICENSE']:
                shutil.copyfile(path,target/path.name)
            elif path.suffix=='.npz':
                keep='returns' not in folder
                inventory.append(dict(folder=folder,file=path.name,bytes=path.stat().st_size,sha256=digest(path),included_in_full=keep))
                if keep:shutil.copyfile(path,target/path.name)
                else:
                    data=np.load(path);mask=data['time']>=data['time'][-1]-50
                    values={key:data[key] for key in data.files if key not in ['time','log_states','physical_states']}
                    values.update(time=data['time'][mask],log_states=data['log_states'][mask],physical_states=data['physical_states'][mask],
                        first_state=data['log_states'][0],final_state=data['log_states'][-1],original_sample_count=len(data['time']))
                    np.savez_compressed(target/(path.stem+'-compact.npz'),**values)
    # Preserve the original broad selection description and add the actual
    # subsequent four-case selection, rather than changing run-time history.
    write(out/'selection-clarifications.json',dict(slow_returns={
        'ids':[415,548,2750,1832],
        'reason':'Selected weak/slow raw-coordinate signals for longer independent return histories; not the strongest raw-rate set.',
        'correction':'The copied generic protocol selection sentence is too broad for this follow-up. It is retained verbatim; this record gives the actual selection.'},
        budget_returns='All four step-budget exits from the260-case long log-coordinate continuation, continued from their saved endpoints.'))
    raw=np.load(args.external/'supramolecular-pilot/complete-or-budgeted.npz')
    log=np.load(args.external/'supramolecular-log-recheck/complete-or-budgeted.npz')
    long=np.load(args.external/'supramolecular-log-long/complete-or-budgeted.npz')
    primary=lambda d,threshold:int(np.sum((d['status']==0)&(d['windows'].min(axis=1)>threshold)))
    budget=json.loads((args.external/'supramolecular-budget-returns/analysis.json').read_text())
    assert len(raw['status'])==len(log['status'])==3394
    assert np.all(log['status']==0) and len(long['status'])==260 and primary(long,.001)==0
    assert budget['completed']==4 and all(r['first_numerical_period_lag']==1 for r in budget['runs'])
    write(out/'audit.json',dict(raw_statuses={str(k):int(np.sum(raw['status']==k)) for k in np.unique(raw['status'])},
        log_statuses={str(k):int(np.sum(log['status']==k)) for k in np.unique(log['status'])},
        long_statuses={str(k):int(np.sum(long['status']==k)) for k in np.unique(long['status'])},
        raw_primary=primary(raw,.005),log_primary=primary(log,.005),log_weak=primary(log,.001),long_retained=primary(long,.001),
        budget_endpoint_numerical_cycles=4,extraction='All solver-root section events and final50tau at original sampling for every return check; full histories hash-indexed locally.',
        conclusion='No retained chaotic-attractor candidate. Finite parameter box and finite observation; no global absence theorem.'))
    write(out/'local-inventory.json',inventory)
    (out/'archive_supramolecular_screen.py').write_bytes(Path(__file__).read_bytes())
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(r['bytes'] for r in entries),files=entries))
    print('files',len(entries),'bytes',sum(r['bytes'] for r in entries),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

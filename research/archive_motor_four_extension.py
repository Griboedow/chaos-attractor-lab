"""Freeze completed long/fresh-start evidence, keeping large raw histories indexed."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);inventory=[]
    folders=['motor-four-100k','motor-four-start-histories','motor-four-start-histories02','motor-four-start-gallery']
    for name in folders:
        source=args.external/name;target=out/name;target.mkdir()
        for path in sorted(source.iterdir()):
            if path.suffix in ['.json','.py','.cl']:shutil.copyfile(path,target/path.name)
            elif path.suffix=='.npz':
                keep=path.name=='display-inputs.npz'
                if keep:shutil.copyfile(path,target/path.name)
                inventory.append(dict(run=name,name=path.name,bytes=path.stat().st_size,sha256=digest(path),included=keep))
    source=args.external/'motor-four-100k/case-3142.npz';d=np.load(source)
    np.savez_compressed(out/'motor-four-100k/sections-and-terminal-segment.npz',
        parameters=d['parameters'],section_times=d['section_times'],section_states=d['section_states'],
        time=d['time'][-12001:],states=d['states'][-12001:],final_state=d['final_state'],final_basis=d['final_basis'],
        log_sums=d['log_sums'],integrated_trace=d['integrated_trace'])
    result=json.loads((out/'motor-four-100k/case-3142.json').read_text())
    assert result['sha256']==digest(source) and result['samples']==1000000 and result['sections']==15543
    assert min(w['spectrum'][0] for w in result['windows'])>.019 and abs(result['spectrum'][1])<2e-5
    reference=json.loads((out/'motor-four-start-histories02/comparisons.json').read_text())
    assert reference['completed']==reference['declared']==32
    for r in reference['runs']:
        assert digest(args.external/'motor-four-start-histories02'/r['file'])==r['sha256']
    a=json.loads((out/'motor-four-start-gallery/analysis.json').read_text())
    assert a['near_permutation_counts']=={'3':32}
    write(out/'local-inventory.json',inventory)
    write(out/'audit.json',dict(long_history_sha256=digest(source),long_samples=1000000,long_sections=15543,
        extracted_history='All exact solver section events plus final1,200tau of uniformly sampled state history; full100ktau history remains hash-indexed locally.',
        fresh_histories=32,completed_comparisons=32,observed_classes=8,near_permutations_per_start=3,
        gallery='Two separately labeled comparison figures; no sum of trajectories presented as an attractor.',
        caveat='No certified attracting set, basin probability, phase locking, topology or historical novelty.'))
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries))
    print(json.dumps(dict(files=len(entries),bytes=sum(e['bytes'] for e in entries),checks='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())

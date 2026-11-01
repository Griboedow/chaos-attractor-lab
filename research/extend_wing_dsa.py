"""Continue all weak wing signals over a longer, predeclared observation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from search_wing_dsa import GPU,write


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=args.source/'complete-or-budgeted.npz';data=np.load(source)
    ids=np.array(json.loads((args.source/'analysis.json').read_text())['weak_ids'])
    p=data['parameters'][ids];initial=data['states'][ids,:4]
    write(out/'protocol.json',dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),source_ids=ids.tolist(),
        selection='Every weak pilot signal, all four windows>.001; independent800tau checks had too few crossings to classify.',
        warm=2000.,observation=40000.,rtol=2e-10,atol=2e-13,maxstep=.03,wall_budget=1200.,step_budget=10000000,
        retention='Every one of four10000tau windows>.001',limitation='Finite continuation, not an asymptotic classification.'))
    np.savez_compressed(out/'declared-inputs.npz',parameters=p,initial=initial,source_ids=ids)
    for name in ['extend_wing_dsa.py','search_wing_dsa.py','wing_dsa.py','wing_dsa.cl']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    result=GPU().scan(p,initial,2000,40000,rtol=2e-10,maxstep=.03);result['source_ids']=ids
    np.savez_compressed(out/'complete-or-budgeted.npz',**result)
    keep=(result['status']==0)&(result['windows'].min(axis=1)>.001)
    write(out/'analysis.json',dict(total=len(ids),statuses={str(k):int(np.sum(result['status']==k)) for k in range(5)},
        retained_source_ids=ids[keep].tolist(),sha256=hashlib.sha256((out/'complete-or-budgeted.npz').read_bytes()).hexdigest(),wall_seconds=result['wall_seconds']))
    print('retained',ids[keep].tolist(),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

"""Long continuation of every weak signal in the positive-coordinate recheck."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from search_rho_adsi import GPU,write
from supramolecular_log import field_jac


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    d=np.load(args.source);ids=np.flatnonzero((d['status']==0)&(d['windows'].min(axis=1)>.001))
    params=d['parameters'][ids];initial=d['states'][ids,:3]
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        source_ids=ids.tolist(),lanes=len(ids),selection='Every completed log-coordinate lane with all4 windows>.001, including all primary>.005 signals.',
        initial='Saved endpoint of the positive-coordinate run; not a fresh start.',warm=500.,observation=8000.,rtol=2e-10,atol=2e-13,maxstep=.02,
        retention='All4 windows>.001; all raw rates retained.',
        limitations='Approximate long-polymer closure; log concentrations may grow without bound or reach physically negligible populations. No novelty or asymptotic-chaos claim.'))
    np.savez_compressed(out/'declared-inputs.npz',source_ids=ids,parameters=params,initial=initial)
    for name in ['extend_supramolecular_log.py','supramolecular_log.py','supramolecular_log.cl','supramolecular.py','search_rho_adsi.py','rho_adsi.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    gpu=GPU(kernel='supramolecular_log.cl',field=field_jac)
    result=gpu.scan(params,initial,500,8000,rtol=2e-10,maxstep=.02,wall=1200);result['source_ids']=ids
    np.savez_compressed(out/'complete-or-budgeted.npz',**result)
    keep=(result['status']==0)&(result['windows'].min(axis=1)>.001)
    write(out/'analysis.json',dict(total=len(ids),statuses={str(i):int((result['status']==i).sum()) for i in range(6)},retained_source_ids=ids[keep].tolist(),
        wall_seconds=result['wall_seconds'],initial_field_error=result['initial_field_error'],
        sha256=hashlib.sha256((out/'complete-or-budgeted.npz').read_bytes()).hexdigest()))
    print('complete',int((result['status']==0).sum()),'retained',ids[keep].tolist(),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

"""Independent close-pair accuracy control for saved motor dispersion starts.

Checks both members, all three perturbation sizes, and four declared phases
per model. This does not assert long chaotic pointwise solver agreement.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from motor_rows import source_force_constraint_field
from motor_rows_gpu import MotorRowsGPU
from audit_motor_permutation_geometry import full


def cpu(job):
    q,p=job
    sol=solve_ivp(lambda t,x:source_force_constraint_field(x,p,4),(0,20),q,
        method='Radau',rtol=2e-12,atol=2e-14,max_step=.01,t_eval=[2.,20.])
    assert sol.success
    return sol.y.T


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    protocol=json.loads((args.source/'protocol.json').read_text())
    data=np.load(args.source/'declared-inputs.npz')
    pair_ids=[i for i,j in enumerate(protocol['jobs']) if j['phase_index'] in [0,11,21,31]]
    ids=np.array([k for i in pair_ids for k in [2*i,2*i+1]])
    q=data['initial'][ids];p=data['parameters'][ids]
    spec=dict(source=str(args.source),source_sha256=hashlib.sha256((args.source/'declared-inputs.npz').read_bytes()).hexdigest(),
        pair_ids=pair_ids,phases=[0,11,21,31],times=[2.,20.],
        independent='Both members integrated with source force-constraint field and Radau, rtol2e-12,atol2e-14,maxstep.01.',
        gpu='Original2e-11/.25 and refined2e-12/.125 settings; each target time integrated afresh from saved inputs.',
        correction='The earlier independent-controls.json tested only one base trajectory per model despite its source comment saying pair. This is the actual paired control.',
        scope='Short-horizon state and separation accuracy, not pointwise agreement after chaotic decorrelation.')
    (out/'protocol.json').write_text(json.dumps(spec,indent=2)+'\n')
    for name in ['check_motor_pair_solver.py','motor_rows.py','motor_rows_gpu.py','motor_rows_gpu.cl','audit_motor_permutation_geometry.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        independent=np.array(list(pool.map(cpu,zip(q,p))))
    gpu=MotorRowsGPU(4);records=[];endpoints=[]
    for setting,(rtol,factor) in enumerate([(2e-11,.25),(2e-12,.125)]):
        for ti,t in enumerate([2.,20.]):
            result=gpu.scan(p,q,transient=0,duration=t,rtol=rtol,step_factor=factor,time_mode='original',max_wall=60)
            assert np.all(result['status']==0)
            states=result['final_states'];endpoints.append(states)
            for k,pair_id in enumerate(pair_ids):
                job=protocol['jobs'][pair_id];source=protocol['sources'][job['source_index']]
                weights=np.array(source['weights']);extent=source['extent']
                a=full(states[2*k:2*k+2],4)*weights/extent
                b=full(independent[2*k:2*k+2,ti],4)*weights/extent
                separation=a[1]-a[0];reference=b[1]-b[0]
                records.append(dict(pair_id=pair_id,**job,time=t,setting=setting,
                    max_state_error=float(np.max(np.linalg.norm(a-b,axis=1))),
                    separation_error=float(np.linalg.norm(separation-reference)),
                    reference_separation=float(np.linalg.norm(reference)),
                    relative_separation_error=float(np.linalg.norm(separation-reference)/np.linalg.norm(reference))))
    np.savez_compressed(out/'controls.npz',initial=q,parameters=p,cpu=independent,gpu=endpoints,pair_ids=pair_ids)
    analysis=dict(records=records,maximum_state_error=max(r['max_state_error'] for r in records),
        maximum_relative_separation_error=max(r['relative_separation_error'] for r in records))
    (out/'analysis.json').write_text(json.dumps(analysis,indent=2)+'\n')
    print({k:v for k,v in analysis.items() if k!='records'},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4);main(p.parse_args())

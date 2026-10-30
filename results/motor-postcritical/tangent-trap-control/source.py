"""Recheck the same captured periodic state with carried and fresh tangents.

Only the diagnostic tangent is changed. The physical initial state and model
are identical. This isolates loss of an expanding tangent component after
floating-point synchronization; it does not repair the underlying trajectory.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

import numpy as np

from motor_rows_gpu import MotorRowsGPU, default_tangent, save, ROOT


def run(output):
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    source=ROOT/'results/motor-postcritical/second-survey/continuation-10000.npz'
    a=np.load(source);lane=25;state=a['final_states'][lane];carried=a['final_tangents'][lane]
    starts=np.repeat(state[None],2,axis=0);p=np.repeat(a['parameters'][lane][None],2,axis=0)
    tangents=np.array([carried,default_tangent(8)])
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source=source.relative_to(ROOT).as_posix(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_lane=lane,state_pair_2_3_differences=[float(state[k]-state[k+1]) for k in [2,10,18]],
        carried_tangent_pair_2_3_differences=[float(carried[k]-carried[k+1]) for k in [2,10,18]],
        original_last_windows=a['lle_windows'][lane].tolist(),
        cases=['same physical state, carried tangent','same physical state, fresh sine tangent'],
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        caveat='A tangent restart is a diagnostic at a captured saddle. It neither changes the physical model nor establishes a chaotic orbit.')
    result=MotorRowsGPU(8).scan(p,starts,tangents=tangents,transient=500.,duration=1000.,
               rtol=1e-10,step_factor=.5,max_steps=8000000,max_wall=180.,time_mode='original')
    result.update(experiment=record)
    save(output/'comparison.npz',result)
    record.update(status='complete' if np.all(result['status']==0) else 'incomplete',
                  lle=result['lle'].tolist(),windows=result['lle_windows'].tolist(),
                  wall_seconds=result['wall_seconds'])
    (output/'summary.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    (output/'source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/motor-postcritical/tangent-trap-control')
    a=parser.parse_args();run(a.output.resolve())

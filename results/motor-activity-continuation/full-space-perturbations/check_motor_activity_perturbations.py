"""Declared small full-space perturbations of selected clustered endpoints."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import numpy as np
import motor_rows
import motor_rows_gpu as gpu_module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose a new output')
    root=Path(__file__).resolve().parents[1]
    source=root/'results/motor-rows/activity-path/N8-lifetime-history.npz'
    old=np.load(source);rng=np.random.default_rng(2026103049)
    parameters=[];states=[];cases=[]
    for sid in [36,42,46,49]:
        index=int(np.flatnonzero(old['source_lane_indices']==sid)[0]);x=old['final_states'][index]
        delta=rng.normal(size=(3,8));delta[2]-=delta[2].mean();delta/=np.linalg.norm(delta)
        reduced=np.r_[delta[0],delta[1],delta[2,:7]]
        for size,sign in [(0.,0),(1e-7,1),(1e-7,-1),(1e-9,1),(1e-9,-1)]:
            y=x+sign*size*reduced
            actual=y-x
            physical_actual=np.r_[actual,-sum(actual[16:])]
            if sign and not np.isclose(np.linalg.norm(physical_actual),size,rtol=1e-5,atol=0):
                raise ArithmeticError('Represented physical perturbation differs from its declared size')
            states.append(y);parameters.append(old['parameters'][index])
            cases.append(dict(source_lane=sid,requested_full_physical_norm=size,sign=sign,
                case='unperturbed state with fresh diagnostic tangent' if not sign else 'declared full-space initial perturbation',
                actual_full_physical_norm=float(np.linalg.norm(physical_actual)),
                exact_reduced_delta=actual.tolist(),physical_unit_direction=delta.tolist()))
    args.output.mkdir(parents=True)
    sources=[Path(__file__),Path(motor_rows.__file__),Path(gpu_module.__file__),Path(gpu_module.__file__).with_suffix('.cl')]
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    for p in sources:(args.output/p.name).write_bytes(p.read_bytes())
    protocol=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_path=source.relative_to(root).as_posix(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_snapshot_sha256=hashes,seed=2026103049,cases=cases,
        selection='Four diagnostic clustered endpoints36,42,46,49, selected before results; no claim of random basin sampling.',
        unchanged_original_cohort=True,limitations=[
            'Initial perturbations are new explicitly declared trajectories, not recovery of erased digits.',
            'Positive tangent growth on an unperturbed saddle is not chaos.',
            'Finite-time departure or persistence does not prove the asymptotic attracting regime.'])
    (args.output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    result=gpu_module.MotorRowsGPU(8).scan(np.array(parameters),np.array(states),transient=1000.,duration=4000.,
        rtol=2e-11,step_factor=.25,max_steps=8000000,max_wall=900.,time_mode='original')
    result.update(source_lane_indices=np.array([c['source_lane'] for c in cases]),protocol=protocol)
    gpu_module.save(args.output/'outcomes.npz',result)


if __name__=='__main__':main()

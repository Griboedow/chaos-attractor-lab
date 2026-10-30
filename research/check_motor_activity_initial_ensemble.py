"""A predeclared finite source-style start ensemble at the clustered candidate."""
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
    if args.output.exists():
        raise FileExistsError('Choose a new output folder')
    root=Path(__file__).resolve().parents[1]
    original_path=root/'results/motor-rows/activity-path/N8.npz'
    original=np.load(original_path)
    parameters=original['parameters'][46]
    seeds=[20251030,20251031,20251032]+list(range(202610300,202610313))
    starts=[];cases=[]
    for i,seed in enumerate(seeds):
        regenerated,positions=motor_rows.source_start(8,seed)
        if i<3:
            initial=original['initial_states'][45+i].copy()
            if not np.array_equal(regenerated,initial):
                raise ArithmeticError('Original source-style seed no longer reproduces exactly')
        else:
            initial=regenerated
        starts.append(initial)
        cases.append(dict(case=i,seed=seed,original_lane=45+i if i<3 else None,
            original_source_initial_state_used=i<3,absolute_initial_positions_nm=positions.tolist()))
    args.output.mkdir(parents=True)
    sources=[Path(__file__),Path(motor_rows.__file__),Path(gpu_module.__file__),Path(gpu_module.__file__).with_suffix('.cl')]
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    for p in sources:
        (args.output/p.name).write_bytes(p.read_bytes())
    protocol=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cases=cases,
        source_path=original_path.relative_to(root).as_posix(),source_sha256=hashlib.sha256(original_path.read_bytes()).hexdigest(),
        source_snapshot_sha256=hashes,parameters=parameters.tolist(),
        prescription='Source equilibrium c=s=1; independently seeded interior absolute shifts uniform[-.01,.01] nm, fixed endpoints, differenced and normalized by2pi/ell.',
        original_cohort_unchanged=True,transient_tau=1000.,observed_tau=9000.,
        selection='All three original tuple15 starts plus13 predeclared additional source-style seeds; no outcome selection.',
        limitations=['Finite targeted ensemble, not a statistically calibrated basin-volume estimate.',
            'Source activity exceeds the illustrated source range; no new experimental calibration.',
            'Positive tangent growth alone cannot distinguish chaotic motion from a periodic saddle.'])
    (args.output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    result=gpu_module.MotorRowsGPU(8).scan(np.tile(parameters,(len(starts),1)),np.array(starts),
        transient=1000.,duration=9000.,rtol=2e-11,step_factor=.25,max_steps=10000000,max_wall=1200.,time_mode='original')
    result['protocol']=protocol
    gpu_module.save(args.output/'outcomes.npz',result)


if __name__=='__main__':
    main()

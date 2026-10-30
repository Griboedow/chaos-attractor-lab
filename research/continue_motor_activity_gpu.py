"""Continue the complete fixed 32-lane activity cohort without reselection."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import motor_rows
import motor_rows_gpu as gpu_module

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--setting',choices=['main','strict'],required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--blocks',type=int,default=40)
    parser.add_argument('--duration',type=float,default=500.)
    parser.add_argument('--wall',type=float,default=2400.)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose a new output directory')
    if args.blocks<1 or min(args.duration,args.wall)<=0:raise ValueError('Positive duration/count/budget required')
    source=ROOT/'results/motor-rows/activity-path/N8-lifetime-history.npz'
    final_block=ROOT/'results/motor-rows/activity-path/N8-lifetime-block-10.npz'
    old=np.load(source);last=np.load(final_block)
    if not np.all(old['status']==0):raise ValueError('The whole preceding cohort must be complete')
    ids=old['source_lane_indices'].copy()
    parameters=old['parameters'].copy()
    states=old['final_states'].copy();tangents=old['final_tangents'].copy()
    if len(ids)!=32 or len(np.unique(ids))!=32:raise ValueError('Expected complete fixed32 cohort')
    for name,value in [('source_lane_indices',ids),('parameters',parameters),('final_states',states),('final_tangents',tangents)]:
        if not np.array_equal(value,last[name]):raise ValueError('Final checkpoint mismatch: '+name)
    rtol,step=(1e-10,.5) if args.setting=='main' else (2e-11,.25)
    args.output.mkdir(parents=True)
    sources=[Path(__file__),Path(motor_rows.__file__),Path(gpu_module.__file__),Path(gpu_module.__file__).with_suffix('.cl')]
    source_hashes={p.name:sha(p) for p in sources}
    for p in sources:(args.output/p.name).write_bytes(p.read_bytes())
    np.savez_compressed(args.output/'initial.npz',source_lane_indices=ids,parameters=parameters,
        states=states,tangents=tangents)
    gpu=gpu_module.MotorRowsGPU(8)
    began=time.perf_counter()
    report=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='running',
        arguments={k:str(v) for k,v in vars(args).items()},source_path=source.relative_to(ROOT).as_posix(),
        source_sha256=sha(source),source_final_block=final_block.relative_to(ROOT).as_posix(),
        source_final_block_sha256=sha(final_block),source_snapshot_sha256=source_hashes,
        exact_saved_state_tangent_parameter_ids_matched=True,cohort_size=32,source_lane_indices=ids.tolist(),
        prior_absolute_tau=7100.,requested_additional_tau=args.blocks*args.duration,
        settings=dict(rtol=rtol,step_factor=step,block_tau=args.duration,window_tau=args.duration/4,
            warmup=0.,max_steps_per_lane=8000000),no_pruning=True,no_reselection=True,
        tangent_policy='Carry the exact stored physical tangent; the existing scanner normalizes its length at each block.',
        device=gpu.device.name,blocks=[],limitations=[
            'Finite-time cohort persistence, not asymptotic chaos or basin probability.',
            'Independent settings start from the same archived endpoints; chaotic trajectories need not remain pointwise close.',
            'Lost positivity is not a certified capture. Incomplete blocks remain unresolved.'])
    history=[];means=[];statuses=[]

    def save():
        report['wall_seconds']=time.perf_counter()-began
        (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        np.savez_compressed(args.output/'history.npz',source_lane_indices=ids,parameters=parameters,
            lle_windows=np.array(history),lle=np.array(means),status=np.array(statuses),
            final_states=states,final_tangents=tangents)

    save()
    for block in range(1,args.blocks+1):
        remaining=args.wall-(time.perf_counter()-began)
        if remaining<=1:
            report['status']='wall budget before next block';break
        initial=states.copy();initial_tangent=tangents.copy()
        result=gpu.scan(parameters,states,tangents=tangents,transient=0.,duration=args.duration,
            rtol=rtol,step_factor=step,max_steps=8000000,max_wall=remaining,time_mode='original')
        result.update(source_lane_indices=ids,carried_initial_tangents=initial_tangent,block_number=block,
            prior_absolute_tau=7100.+(block-1)*args.duration,no_pruning=True,no_warmup=True)
        gpu_module.save(args.output/f'block-{block:02d}.npz',result)
        valid=result['status']==0
        positive=valid&(result['lle']>.005)&np.all(result['lle_windows']>.001,axis=1)
        entry=dict(block=block,completed=int(valid.sum()),positive_screen=int(positive.sum()),
            positive_source_lanes=ids[positive].tolist(),nonpositive_complete_source_lanes=ids[valid&~positive].tolist(),
            status_counts={str(i):int(np.sum(result['status']==i)) for i in range(5)},
            exact_initial_state_copy=bool(np.array_equal(initial,result['initial_states'])),
            wall_seconds=result['wall_seconds'],common_additional_completed_tau=(block if valid.all() else block-1)*args.duration)
        history.append(result['lle_windows'].copy());means.append(result['lle'].copy());statuses.append(result['status'].copy())
        report['blocks'].append(entry)
        states=result['final_states'].copy();tangents=result['final_tangents'].copy()
        report['status']='complete' if block==args.blocks and valid.all() else 'checkpoint'
        save();print(json.dumps(entry),flush=True)
        if not valid.all():
            report['status']='incomplete block preserved; cohort not pruned';break
    if report['status']=='complete':
        retained=np.all(np.array(statuses)==0,axis=0)&np.all(np.array(means)>.005,axis=0)&np.all(np.array(history)>.001,axis=(0,2))
        report['retained_all_new_windows_source_lanes']=ids[retained].tolist()
        report['retained_all_new_windows_count']=int(retained.sum())
    save()


if __name__=='__main__':main()

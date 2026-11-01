"""Declared one-parameter slice to investigate the motor braid's mechanism."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from motor_rows_gpu import MotorRowsGPU,save


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);anchor=source['parameters'];assert source['states'].shape[1]==11
    b=np.linspace(.8,2.2,281);tuples=np.column_stack([np.full(len(b),anchor[0]),np.full(len(b),anchor[1]),b])
    parameters=np.repeat(tuples,2,axis=0);initial=np.zeros((len(parameters),11))
    initial[:,:4]=1;initial[:,4:8]=parameters[:,2,None]
    starts=np.random.default_rng(2026110101).uniform(-.05,.05,(2,4));starts-=starts.mean(axis=1)[:,None]
    initial[:,8:]=np.tile(starts[:,:3],(len(tuples),1))
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        N=4,dimension=11,anchor_parameters=anchor.tolist(),b=b.tolist(),lanes=len(parameters),seed=2026110101,
        initial='Two declared fresh zero-sum displacement vectors, c=1 and s=b for each b. Not numerical parameter continuation.',
        warm=1000.,observation=4000.,rtol=1e-9,time='original tau',
        scope='One-parameter slice for mechanism/periodic-window follow-up. No new-model, topology or novelty claim.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf8')
    np.savez_compressed(out/'declared-inputs.npz',parameters=parameters,initial=initial,tuple_id=np.repeat(np.arange(len(tuples)),2),start_id=np.tile([0,1],len(tuples)))
    for name in ['search_motor_four_b_slice.py','motor_rows_gpu.py','motor_rows_gpu.cl','motor_rows.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    gpu=MotorRowsGPU(4)
    result=gpu.scan(parameters,initial,transient=1000,duration=4000,rtol=1e-9,time_mode='original',step_factor=.5,max_steps=15000000,max_wall=1200)
    result.update(tuple_id=np.repeat(np.arange(len(tuples)),2),start_id=np.tile([0,1],len(tuples)),
        protocol_sha256=hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest())
    save(out/'complete-or-budgeted.npz',result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

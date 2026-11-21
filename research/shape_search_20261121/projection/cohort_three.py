"""Estimate three-rate persistence across 1,024 independent starts per retained map."""
import argparse
import json
from pathlib import Path
import numpy as np
import pyopencl as cl

from third_growth import kernel_for,run_gpu
from model import replay
from reproduce import raw_parameters


def three(windows):
    return bool(np.isfinite(windows).all() and windows.min()>.001)


def main(args):
    args.cache.mkdir(parents=True,exist_ok=True)
    dev=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU)
    ctx=cl.Context([dev]);q=cl.CommandQueue(ctx);kernel=kernel_for(ctx,'single')
    rows=[]
    for path in sorted(args.source.glob('case-*.json')):
        rec=json.loads(path.read_text())
        if not rec['three_positive_passes']:continue
        raw=raw_parameters(rec);case=rec['case']
        rng=np.random.default_rng(84421+case)
        states=rng.uniform(-2,2,(args.starts,3)).astype('f4')
        checks=rng.choice(args.starts,size=args.checks,replace=False)
        checkpoint=args.cache/f'cohort-{case}.npz'
        if checkpoint.exists():
            with np.load(checkpoint) as f:
                assert np.array_equal(f['initial'],states)
                gpu=f['gpu'];reeval=f['cpu']
        else:
            data=np.broadcast_to(raw[:,None,:],(len(raw),args.starts,4)).astype('f4').copy()
            data[1,:,:3]=states
            gpu=run_gpu(q,kernel,data,'single',32768,8192)
            reevaluation=[]
            for i in checks:
                _,windows,_,_,_=replay(raw,states[i].astype(float),65536,32768,4,0,True)
                reevaluation.append(windows)
            reeval=np.array(reevaluation)
            np.savez_compressed(checkpoint,initial=states,gpu=gpu,check_indices=checks,cpu=reeval)
        flagged=np.array([three(v) for v in gpu]);truth=np.array([three(v) for v in reeval])
        counts=np.bincount(2*flagged[checks].astype(int)+truth.astype(int),minlength=4)
        rows.append(dict(case=case,seed=84421+case,starts=args.starts,gpu_three=int(flagged.sum()),cpu_checked=args.checks,cpu_three=int(truth.sum()),paired_counts=counts.tolist()))
        print(case,int(flagged.sum()),int(truth.sum()),flush=True)
    doc=dict(threshold='All three exponents above 0.001 in each of four windows',volume='Uniform random initial states in [-2,2]^3; independently seeded for each case; identical binary32 input states for GPU and widened binary64 CPU.',gpu_burn=32768,gpu_window=8192,cpu_burn=65536,cpu_window=32768,paired_order=['neither','CPU only','GPU only','both'],checkpoint_prefix='projection/three-cohort/',rows=rows)
    args.out.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cache',type=Path,required=True);p.add_argument('--starts',type=int,default=1024);p.add_argument('--checks',type=int,default=32)
    main(p.parse_args())

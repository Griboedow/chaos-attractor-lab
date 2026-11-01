"""Continue all forty flagged pilot cases, including those already decaying."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import set_num_threads
import model
from followup import follow

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--block',type=int,default=250000);ap.add_argument('--blocks',type=int,default=20)
    ap.add_argument('--keep',type=int,default=25000)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    for n in (3,4):
        parent=args.input/f'n{n}.json';data=json.loads(parent.read_text());rows=data['runs']
        L=np.array([r['log_kernel'] for r in rows]);z=np.array([r['blocks'][-1]['checkpoint'] for r in rows]);w=np.array([r['omega'] for r in rows])
        t=time.perf_counter();stats,checkpoints,tails=follow(L,z,w,args.block,args.blocks,args.keep)
        output=[]
        for i,r in enumerate(rows):
            start=r['blocks'][-1]['stop_step']
            output.append(dict(case=r['case'],shape=[n,n],omega=r['omega'],start_step=start,
                starting_state=z[i].tolist(),log_kernel=r['log_kernel'],
                blocks=[dict(stop_step=start+(j+1)*args.block,rate=float(s[0]),min_abs_state=float(s[1]),max_abs_state=float(s[2]),
                    final_step_residual=float(s[3]),mean=s[4:4+2*n-2].tolist(),std=s[4+2*n-2:].tolist(),checkpoint=checkpoints[i,j].tolist()) for j,s in enumerate(stats[i])]))
        result=dict(block_length=args.block,blocks=args.blocks,retained_final_window_steps=args.keep,
            elapsed_seconds=time.perf_counter()-t,parent_sha256=hashlib.sha256(parent.read_bytes()).hexdigest(),runs=output,
            sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(model.__file__),Path('followup.py'))})
        (args.output/f'n{n}.json').write_text(json.dumps(result,indent=2)+'\n')
        np.savez_compressed(args.output/f'n{n}-final-{args.keep}-steps.npz',case=np.array([r['case'] for r in rows]),states=tails)
        print(json.dumps(dict(size=n,elapsed=result['elapsed_seconds'],cases=[dict(case=r['case'],first=r['blocks'][0]['rate'],last=r['blocks'][-1]['rate'],norm=r['blocks'][-1]['max_abs_state']) for r in output])),flush=True)

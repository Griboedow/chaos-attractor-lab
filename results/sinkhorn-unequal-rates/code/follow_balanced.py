"""Continue every balanced-kernel flag with the actual pilot origin."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import set_num_threads
from follow_dual import follow

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--continuation',action='store_true');ap.add_argument('--block',type=int,default=50000);ap.add_argument('--blocks',type=int,default=20)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    parent=args.input/'3x2.json';data=json.loads(parent.read_text())
    if args.continuation:
        rows=[r for r in data['runs'] if r['completed_steps']==data['block_length']*data['blocks']];starts=np.array([r['blocks'][-1]['checkpoint'] for r in rows]);origin=data['stop_step']
    else:
        rows=[r for r in data['runs'] if r['rates'] and r['rates']['mean']>.005];origin=data['warm']+data['observed']
        ends=np.load(args.input/'3x2-endpoints-and-last2048.npz')['endpoints'];starts=np.array([ends[r['case']] for r in rows])
    L=np.array([r['log_kernel'] for r in rows]);ws=np.array([[r['wx'],r['wy']] for r in rows]);t=time.perf_counter()
    stats,checkpoints,tails,completed=follow(L,starts,ws,args.block,args.blocks,50000);out=[]
    for i,r in enumerate(rows):
        out.append(dict(case=r['case'],shape=[3,2],wx=r['wx'],wy=r['wy'],log_kernel=r['log_kernel'],starting_state=starts[i].tolist(),completed_steps=int(completed[i]),
            blocks=[dict(stop_step=origin+(j+1)*args.block,rate=float(s[0]),min_abs_state=float(s[1]),max_abs_state=float(s[2]),final_step_residual=float(s[3]),mean=s[4:7].tolist(),std=s[7:].tolist(),checkpoint=checkpoints[i,j].tolist()) for j,s in enumerate(stats[i]) if np.isfinite(s).all()]))
    result=dict(runs=out,start_step=origin,stop_step=origin+args.block*args.blocks,block_length=args.block,blocks=args.blocks,retained_final_window_steps=50000,parent_sha256=hashlib.sha256(parent.read_bytes()).hexdigest(),elapsed_seconds=time.perf_counter()-t,
        sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('model.py'),Path('model_dual.py'),Path('follow_dual.py'))})
    (args.output/'3x2.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');np.savez_compressed(args.output/'3x2-final-50000-states.npz',case=np.array([r['case'] for r in rows]),states=tails)
    print(json.dumps(dict(flags=len(rows),complete=int(sum(completed==args.block*args.blocks)),elapsed=result['elapsed_seconds'],cases=[dict(case=r['case'],completed=r['completed_steps'],last=r['blocks'][-1]['rate'] if r['blocks'] else None) for r in out])),flush=True)

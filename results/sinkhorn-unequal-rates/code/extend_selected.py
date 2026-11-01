"""A further 100 million steps for the two selected balanced-kernel cases."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import set_num_threads
from follow_dual import follow

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(2);data=json.loads(args.input.read_text());rows=[r for r in data['runs'] if r['case'] in (150,987)];assert len(rows)==2
    origin=data['stop_step'];starts=np.array([r['blocks'][-1]['checkpoint'] for r in rows]);L=np.array([r['log_kernel'] for r in rows]);ws=np.array([[r['wx'],r['wy']] for r in rows]);t=time.perf_counter()
    stats,checkpoints,tails,completed=follow(L,starts,ws,5000000,20,250000);out=[]
    for i,r in enumerate(rows):
        out.append(dict(case=r['case'],shape=[3,2],wx=r['wx'],wy=r['wy'],log_kernel=r['log_kernel'],starting_state=starts[i].tolist(),completed_steps=int(completed[i]),blocks=[dict(stop_step=origin+(j+1)*5000000,rate=float(s[0]),min_abs_state=float(s[1]),max_abs_state=float(s[2]),final_step_residual=float(s[3]),mean=s[4:7].tolist(),std=s[7:].tolist(),checkpoint=checkpoints[i,j].tolist()) for j,s in enumerate(stats[i]) if np.isfinite(s).all()]))
    result=dict(runs=out,start_step=origin,stop_step=origin+100000000,block_length=5000000,blocks=20,retained_final_window_steps=250000,parent_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),elapsed_seconds=time.perf_counter()-t,sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('model.py'),Path('model_dual.py'),Path('follow_dual.py'))})
    (args.output/'3x2.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');np.savez_compressed(args.output/'3x2-final-250000-states.npz',case=np.array([r['case'] for r in rows]),states=tails)
    print(json.dumps(dict(elapsed=result['elapsed_seconds'],cases=[dict(case=r['case'],completed=r['completed_steps'],last=r['blocks'][-1]['rate'] if r['blocks'] else None) for r in out])),flush=True)

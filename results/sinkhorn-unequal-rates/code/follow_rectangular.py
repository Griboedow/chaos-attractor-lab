"""Continue every completed-window rectangular pilot flag."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import set_num_threads
from base_followup import follow

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--pilot',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    for m,n in [(2,3),(3,2),(2,4),(4,2)]:
        parent=args.pilot/f'{m}x{n}.json';pilot=json.loads(parent.read_text());rows=[r for r in pilot['runs'] if not r['capture_step'] and r['finite_rate']>1e-5]
        ends=np.load(args.pilot/f'{m}x{n}-endpoints-and-last2048.npz')['endpoints'];starts=np.array([ends[r['case']] for r in rows]);L=np.array([r['log_kernel'] for r in rows]);w=np.array([r['omega'] for r in rows]);t=time.perf_counter()
        stats,checkpoints,tails=follow(L,starts,w,50000,20,50000);d=m+n-2;out=[]
        for i,r in enumerate(rows):
            out.append(dict(case=r['case'],shape=[m,n],omega=r['omega'],log_kernel=r['log_kernel'],start_step=60000,starting_state=starts[i].tolist(),
                blocks=[dict(stop_step=60000+(j+1)*50000,rate=float(s[0]),min_abs_state=float(s[1]),max_abs_state=float(s[2]),final_step_residual=float(s[3]),mean=s[4:4+d].tolist(),std=s[4+d:].tolist(),checkpoint=checkpoints[i,j].tolist()) for j,s in enumerate(stats[i])]))
        result=dict(runs=out,block_length=50000,blocks=20,retained_final_window_steps=50000,parent_sha256=hashlib.sha256(parent.read_bytes()).hexdigest(),elapsed_seconds=time.perf_counter()-t,
            sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('model.py'),Path('base_followup.py'))})
        (args.output/f'{m}x{n}.json').write_text(json.dumps(result,indent=2)+'\n')
        np.savez_compressed(args.output/f'{m}x{n}-final-50000-states.npz',case=np.array([r['case'] for r in rows]),states=tails)
        print(json.dumps(dict(shape=[m,n],elapsed=result['elapsed_seconds'],cases=[dict(case=r['case'],first=r['blocks'][0]['rate'],last=r['blocks'][-1]['rate'],norm=r['blocks'][-1]['max_abs_state']) for r in out])),flush=True)

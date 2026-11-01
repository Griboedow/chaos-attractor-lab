"""Unmodified installed SciPy RK23 trajectories, including its actual clock."""
import argparse,hashlib,inspect,json,time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import RK23
import scipy.integrate._ivp.rk as upstream
from search_rk23_controller import rk23_step

def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    for n in ['check_rk23_upstream.py','search_rk23_controller.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    records=[];source=Path(upstream.__file__)
    for ident in args.ids:
        angle,rtol=d['params'][ident];lam=-np.cos(angle)+1j*np.sin(angle)
        solver=RK23(lambda t,y:lam*y,0.,np.array([1.+0j]),1e100,rtol=rtol,atol=1.,first_step=.1)
        states=np.empty((400001,5));states[0]=[1.,0.,.1,0.,0.]
        start=time.time();first=[]
        for n in range(400000):
            old=np.array([np.log(solver.h_abs),np.log(abs(solver.y[0]))])
            if n<20:expected=rk23_step(old,np.array([angle,rtol]))[0]
            solver.step();assert solver.status=='running'
            z=solver.y[0];states[n+1]=[z.real,z.imag,solver.h_abs,solver.t,solver.h_previous]
            if n<20:first.append(float(np.max(abs(np.array([np.log(solver.h_abs),np.log(abs(z))])-expected))))
        np.savez_compressed(args.out/f'case-{ident}.npz',states=states,parameters=[angle,rtol])
        record=dict(id=ident,parameters=[angle,rtol],atol=1.,first_step=.1,steps=400000,seconds=time.time()-start,
                    first20_max_one_step_log_residual=max(first),evaluations=solver.nfev,
                    columns=['Re y','Im y','proposed h','time','previous accepted h'],scipy_version=scipy.__version__,
                    installed_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
        records.append(record);(args.out/'report.json').write_text(json.dumps(dict(runs=records),indent=2),encoding='utf-8')
        print(json.dumps(record),flush=True)
    (args.out/'installed-rk.py').write_bytes(source.read_bytes())
    for n in ['check_rk23_upstream.py','search_rk23_controller.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',default=[11,396,719,1024]);main(p.parse_args())

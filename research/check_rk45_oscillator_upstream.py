"""Run unmodified SciPy on a physical initial displacement, then audit its map."""
import argparse,hashlib,json,time,shutil
from pathlib import Path
import numpy as np
from numba import njit
from scipy.integrate import RK45
import scipy,scipy.integrate._ivp.rk as upstream
from search_rk45_oscillator import matrices,oscillator_step
from follow_rk45_oscillator import orthogonalize

@njit(cache=True)
def analyze(states,params):
    ms=matrices(params[0]);q,_,_=orthogonalize(np.array([[1.,.2,-.3],[.4,1.,.15],[-.25,.35,1.]]))
    rates=np.zeros((29,3));residual=0.;detlogs=np.zeros(29);rejections=0
    for n in range(100000,400000):
        s=np.array([states[n,0],states[n,1],np.log(states[n,2])])
        expected,j,reject,status=oscillator_step(s,params[1],*ms)
        actual=np.array([states[n+1,0],states[n+1,1],np.log(states[n+1,2])])
        if status:raise ValueError('Transcription failed')
        residual=max(residual,np.max(np.abs(expected-actual)))
        q,r,status=orthogonalize(j@q)
        if status:raise ValueError('QR failed')
        if n>=110000:
            k=(n-110000)//10000;rates[k]+=r;rejections+=reject
            det=j[0,0]*(j[1,1]*j[2,2]-j[1,2]*j[2,1])-j[0,1]*(j[1,0]*j[2,2]-j[1,2]*j[2,0])+j[0,2]*(j[1,0]*j[2,1]-j[1,1]*j[2,0])
            detlogs[k]+=np.log(abs(det))
    return rates/10000,detlogs/10000,residual,rejections

def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    for n in ['check_rk45_oscillator_upstream.py','search_rk45_oscillator.py','follow_rk45_oscillator.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    path=Path(upstream.__file__);(args.out/'installed-rk.py').write_bytes(path.read_bytes())
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    records=[]
    for ident in args.ids:
        params=d['params'][ident];A=np.array([[0.,1.],[-1.,-2*params[0]]])
        begin=time.time();reused=args.reuse/f'case-{ident}.npz' if args.reuse else None
        evaluations=None
        if reused is not None and reused.exists():
            with np.load(reused) as saved:
                assert np.array_equal(saved['params'],params) and float(saved['atol'])==args.atol
                states=saved['states']
            assert states.shape==(400001,5) and np.all(np.isfinite(states))
            shutil.copy2(reused,args.out/reused.name)
        else:
            solver=RK45(lambda t,y:A@y,0.,np.array([1.,0.]),1e100,rtol=params[1],atol=args.atol,first_step=.1)
            states=np.empty((400001,5));states[0]=[1/args.atol,0.,.1,0.,0.]
            for n in range(400000):
                solver.step();assert solver.status=='running'
                states[n+1]=[solver.y[0]/args.atol,solver.y[1]/args.atol,solver.h_abs,solver.t,solver.h_previous]
            evaluations=solver.nfev
            np.savez_compressed(args.out/f'case-{ident}.npz',states=states,params=params,atol=args.atol)
        rates,detlogs,residual,reject=analyze(states,params)
        assert residual<1e-6,residual
        np.savez_compressed(args.out/f'case-{ident}-analysis.npz',rates=rates,detlogs=detlogs)
        r=dict(id=ident,params=params.tolist(),atol=args.atol,physical_initial=[1.,0.],first_step=.1,
               steps=400000,warm=100000,alignment=10000,observed=290000,windows=29,means=rates.mean(0).tolist(),minima=rates.min(0).tolist(),
               max_one_step_residual=residual,max_window_determinant_identity_residual=float(np.max(abs(rates.sum(1)-detlogs))),
               predicted_late_rejections=int(reject),evaluations=evaluations,seconds_this_run=time.time()-begin,
               reused_trajectory=str(reused) if reused is not None and reused.exists() else None,
               scipy_version=scipy.__version__,installed_source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
               columns=['position / atol','velocity / atol','proposed step','time','previous accepted step'],
               scope='Unmodified installed solver, including its floating-point clock. Tangent rates use the mathematical transcription evaluated on these actual states.')
        records.append(r);(args.out/'report.json').write_text(json.dumps(dict(runs=records),indent=2),encoding='utf-8');print(json.dumps(r),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--reuse',type=Path);p.add_argument('--ids',type=int,nargs='+',default=[1890,335,448,2006]);p.add_argument('--atol',type=float,default=1e-6);main(p.parse_args())

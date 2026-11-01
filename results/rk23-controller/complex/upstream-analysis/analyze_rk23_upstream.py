"""Audit saved actual SciPy complex-scalar trajectories, without rerunning them."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit
from search_rk23_controller import rk23_step

@njit(cache=True)
def analyze(states,params):
    q=np.array([[.9,-.435889894],[.435889894,.9]]);q/=np.sqrt(np.sum(q[:,0]**2))
    rates=np.zeros((29,2));residual=0.;detlogs=np.zeros(29)
    for n in range(100000,400000):
        s=np.array([np.log(states[n,2]),np.log(np.hypot(states[n,0],states[n,1]))])
        expected,j,_,_,status=rk23_step(s,params)
        if status:raise ValueError('Chart or map failure')
        actual=np.array([np.log(states[n+1,2]),np.log(np.hypot(states[n+1,0],states[n+1,1]))])
        residual=max(residual,np.max(np.abs(expected-actual)))
        jq=j@q;v=jq[:,0];r0=np.sqrt(np.sum(v*v));q[:,0]=v/r0
        v=jq[:,1]-np.dot(q[:,0],jq[:,1])*q[:,0];v-=np.dot(q[:,0],v)*q[:,0];r1=np.sqrt(np.sum(v*v));q[:,1]=v/r1
        if n>=110000:
            k=(n-110000)//10000;rates[k,0]+=np.log(r0);rates[k,1]+=np.log(r1)
            detlogs[k]+=np.log(abs(j[0,0]*j[1,1]-j[0,1]*j[1,0]))
    return rates/10000,detlogs/10000,residual

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);rows=[]
    for ident in [11,396,719,1024]:
        source=(args.control if ident==1024 else args.source)/f'case-{ident}.npz'
        with np.load(source) as f:states=f['states'];params=f['parameters']
        rates,detlogs,residual=analyze(states,params);assert residual<1e-6
        np.savez_compressed(args.out/f'case-{ident}.npz',rates=rates,detlogs=detlogs)
        row=dict(id=ident,params=params.tolist(),source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),warm=100000,alignment=10000,observed=290000,windows=29,means=rates.mean(0).tolist(),minima=rates.min(0).tolist(),max_one_step_log_residual=residual,max_determinant_identity_residual=float(np.max(np.abs(rates.sum(1)-detlogs))))
        rows.append(row);print(json.dumps(row),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(runs=rows,scope='Two-dimensional amplitude/step quotient only. Exact phase symmetry adds a neutral direction; these are not three-dimensional novelty claims.'),indent=2),encoding='utf-8')
    for n in ['analyze_rk23_upstream.py','search_rk23_controller.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--control',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

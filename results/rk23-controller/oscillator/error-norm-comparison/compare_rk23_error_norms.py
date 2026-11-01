"""Same oscillator eigenvalues, a deliberately different rotation-invariant error rule."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit
from search_rk23_controller import rk23_step

@njit(cache=True)
def run(params,amplitude):
    s=np.array([np.log(.1),np.log(amplitude)]);v=np.array([.9,.435889894]);v/=np.sqrt(np.sum(v*v))
    rates=np.zeros((20,2));tail=np.empty((4096,2));bounds=np.array([[np.inf,np.inf],[-np.inf,-np.inf]])
    rejects=0
    for n in range(1100000):
        s,J,_,r,status=rk23_step(s,params)
        if status or not np.all(np.isfinite(s)) or np.max(np.abs(s))>100:raise ValueError('Map chart failure')
        v=J@v;norm=np.sqrt(np.sum(v*v));v/=norm;lead=np.log(norm);det=J[0,0]*J[1,1]-J[0,1]*J[1,0]
        if n>=100000:
            k=(n-100000)//50000;rates[k,0]+=lead;rates[k,1]+=np.log(abs(det))-lead;rejects+=r
            bounds[0]=np.minimum(bounds[0],s);bounds[1]=np.maximum(bounds[1],s)
        if n>=1100000-4096:tail[n-(1100000-4096)]=s
    return rates/50000,tail,bounds,rejects

def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    with np.load(args.source) as f:params=f['params']
    rows=[]
    for ident in [381,354,883,37,526,689]:
        zeta,rtol=params[ident];omega=np.sqrt(1-zeta*zeta);p=np.array([np.arccos(zeta),rtol]);amplitude=1/omega;begin=time.time()
        rates,tail,bounds,rejects=run(p,amplitude);period=0
        for k in range(1,513):
            if np.max(np.abs(tail[k:]-tail[:-k]))<1e-9:period=k;break
        np.savez_compressed(args.out/f'case-{ident}.npz',params=p,oscillator_params=params[ident],rates=rates,tail=tail,bounds=bounds)
        row=dict(id=ident,oscillator_params=params[ident].tolist(),complex_params=p.tolist(),complex_initial_amplitude=amplitude,means=rates.mean(0).tolist(),minima=rates.min(0).tolist(),quotient_period=period,rejections=int(rejects),seconds=time.time()-begin)
        rows.append(row);print(json.dumps(row),flush=True)
    report=dict(runs=rows,source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),warm=100000,observed=1000000,windows=20,
                transformation='z=q-i*(p+zeta*q)/sqrt(1-zeta^2); then z_dot=(-zeta+i*sqrt(1-zeta^2))*z. Initial physical (q,p)=(1,0) has |z|=1/sqrt(1-zeta^2).',
                changed_rule='Uses |embedded complex error| / [1+rtol*max(|z_old|,|z_new|)], rather than original-axis componentwise scales and RMS.',
                scope='Error-norm ablation, intentionally a different numerical algorithm for a conjugate linear ODE. Exact phase symmetry contributes a neutral direction when amplitudes remain nonzero and bounded. Not an invariant-set inequivalence proof or a new coordinate-dependence discovery.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    for n in ['compare_rk23_error_norms.py','search_rk23_controller.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

"""Full tangents along saved, actual TensorFlow Addons histories."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit
from yogi_tanh_map import step_tangent


@njit
def along(states,p):
    q,_=np.linalg.qr(np.array([[1.,2.,3.],[4.,5.,7.],[2.,6.,8.]]))
    acc=np.zeros(3);dacc=0.;windows=np.zeros((19,3));detwindows=np.zeros(19);error=0.
    for t in range(199999):
        expected,q,det=step_tangent(states[t],q,p)
        error=max(error,np.max(abs(expected-states[t+1])/(1+abs(states[t+1]))))
        logs=np.empty(3)
        for j in range(3):
            for _ in range(2):
                for h in range(j):q[:,j]-=np.dot(q[:,j],q[:,h])*q[:,h]
            norm=np.linalg.norm(q[:,j]);q[:,j]/=norm;logs[j]=np.log(norm)
        # 9,999 transitions for alignment, then 19 complete10k windows.
        if t>=9999:
            acc+=logs;dacc+=np.log(abs(det));n=t-9999+1
            if n%10000==0:windows[n//10000-1]=acc/10000;detwindows[n//10000-1]=dacc/10000;acc[:]=0.;dacc=0.
    return windows,detwindows,error


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);a=json.loads((args.source/'analysis.json').read_text());rows=[];allw=[];alld=[];ids=[]
    for original in a['runs']:
        ident=original['id'];file=args.source/f'case-{ident}.npz'
        with np.load(file) as f:d={k:f[k] for k in f.files}
        states=d['states'];assert states.shape==(300000,3)
        w,det,err=along(states[-200000:],d['parameters']);assert np.isfinite(w).all() and err<1e-10
        residual=float(np.max(abs(w.sum(1)-det)));assert residual<1e-7
        tail=states[-4096:];tol=1e-9*(1+np.max(abs(tail)));period=next((lag for lag in range(1,129) if np.max(abs(tail[lag:]-tail[:-lag]))<tol),None)
        fixed_rates=None
        if np.max(abs(states[-1,:2]))<1e-100:
            b,c,e,k=d['parameters'];den=np.sqrt(states[-1,2])+e
            J=np.array([[1-(1-b)/den,-b/den,0.],[1-b,b,0.],[0.,0.,1.]])
            fixed_rates=sorted(np.log(abs(np.linalg.eigvals(J))).tolist(),reverse=True)
            assert np.max(abs(np.sort(w.mean(0))[::-1]-fixed_rates))<1e-4
        r=dict(id=ident,parameters=d['parameters'].tolist(),spectrum=w.mean(0).tolist(),window_minima=w.min(0).tolist(),window_maxima=w.max(0).tolist(),
            period_up_to128=period,exact_fixed_line_linear_rates=fixed_rates,step_max_scaled_residual=float(err),determinant_sum_max_residual=residual,source_sha256=hashlib.sha256(file.read_bytes()).hexdigest())
        rows.append(r);ids.append(ident);allw.append(w);alld.append(det);print(r,flush=True)
    np.savez_compressed(out/'windows.npz',ids=ids,windows=allw,determinant_windows=alld)
    (out/'analysis.json').write_text(json.dumps(dict(source=str(args.source),runs=rows,completed=len(rows),declared=len(a['runs']),warmup=100000,alignment_transitions=9999,observed_transitions=190000,
        scope='Actual library states with analytic full tangents at the exact effective library coefficients. No inference of equality with a different autonomous basin.'),indent=2)+'\n')
    for n in ['analyze_yogi_tanh_upstream.py','yogi_tanh_map.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

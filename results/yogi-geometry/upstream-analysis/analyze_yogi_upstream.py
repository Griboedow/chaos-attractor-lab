"""Evaluate the normalized recurrence and full branch tangent along actual Optax states."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
from yogi_map import step_tangent


@njit
def along_actual(states,p):
    q=np.eye(3);logs=np.zeros(3);windows=np.empty((19,3));error=0.;margin=np.inf
    for t in range(200000):
        expected,q=step_tangent(states[t],q,p)
        error=max(error,np.max(np.abs(expected-states[t+1])/(1+np.abs(states[t+1]))))
        growth=np.empty(3)
        for j in range(3):
            for _ in range(2):
                for k in range(j):q[:,j]-=np.dot(q[:,j],q[:,k])*q[:,k]
            norm=np.linalg.norm(q[:,j]);q[:,j]/=norm;growth[j]=np.log(norm)
        if t>=10000:
            logs+=growth;margin=min(margin,abs(states[t,2]-states[t,0]**2)/(1+states[t,2]))
            n=t-10000+1
            if n%10000==0:windows[n//10000-1]=logs/10000;logs[:]=0.
    return windows,error,margin


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);reference=np.load(args.reference);rows=[];all_windows=[]
    for i,ident in enumerate(source['ids']):
        states=np.vstack([source['warm_endpoint'][i],source['states'][:,i]]);p=source['parameters'][i]
        windows,error,margin=along_actual(states,p);assert error<1e-12 and np.isfinite(windows).all()
        sum_error=float(np.max(abs(windows.sum(1)-np.log(p[0]))));assert sum_error<1e-10
        all_windows.append(windows);tail=states[-4096:]
        returns=[float(np.max(np.linalg.norm(tail[lag:]-tail[:-lag],axis=1))) for lag in range(1,129)]
        period=next((k+1 for k,r in enumerate(returns) if r<1e-9*(1+np.max(abs(tail)))),None)
        k=list(reference['ids']).index(ident);r=reference['tails'][k];extent=np.linalg.norm(np.ptp(r,axis=0))
        unique,inverse=np.unique(states[-20000:],axis=0,return_inverse=True)
        distances=cKDTree(r).query(unique)[0][inverse]/max(extent,1e-12)
        row=dict(id=int(ident),parameters=p.tolist(),autonomous_step_max_scaled_residual=float(error),spectrum=windows.mean(0).tolist(),
            window_ranges=np.array([windows.min(0),windows.max(0)]).tolist(),windows=windows.tolist(),determinant_sum_max_residual=sum_error,
            observed_after_alignment=190000,minimum_scaled_switch_margin=float(margin),tail_period_up_to128=period,
            tail_extents=np.ptp(tail,axis=0).tolist(),final_state=states[-1].tolist(),
            support_distances_q50_q95_q99_max=np.quantile(distances,[.5,.95,.99,1]).tolist(),support_normalization=max(extent,1e-12))
        rows.append(row);print(ident,row['spectrum'],'period',period,'map residual',error,flush=True)
    np.savez_compressed(out/'windows.npz',ids=source['ids'],windows=all_windows)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows,completed=len(rows),declared=len(source['ids']),
        source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        reference=str(args.reference),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        scope='All states come from actual Optax0.2.8. Analytic branch tangents are evaluated along them. Independent scalar finite-separation checks are separate.'),indent=2)+'\n')
    for n in ['analyze_yogi_upstream.py','yogi_map.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

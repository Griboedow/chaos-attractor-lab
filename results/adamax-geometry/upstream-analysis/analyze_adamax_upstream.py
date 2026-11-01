"""Check actual upstream state updates and tangent growth along those states."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
from adamax_map import step_tangent


@njit
def along_actual(states,p):
    q=np.eye(3)[:,:2].copy();logs=np.zeros(2);windows=np.empty((19,2));error=0.;resets=0
    for t in range(200000):
        expected,q,reset=step_tangent(states[t],q,p)
        error=max(error,np.max(np.abs(expected-states[t+1])/(1+np.abs(states[t+1]))))
        growth=np.empty(2)
        for j in range(2):
            for _ in range(2):
                for k in range(j):q[:,j]-=np.dot(q[:,j],q[:,k])*q[:,k]
            norm=np.linalg.norm(q[:,j]);q[:,j]/=norm;growth[j]=np.log(norm)
        if t>=10000:
            logs+=growth;resets+=int(reset)
            n=t-10000+1
            if n%10000==0:windows[n//10000-1]=logs/10000;logs[:]=0.
    return windows,error,resets


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);reference=np.load(args.reference);rows=[];all_windows=[]
    for i,ident in enumerate(source['ids']):
        states=source['states'][i];p=source['parameters'][i]
        windows,error,resets=along_actual(states,p);assert error<1e-12 and np.isfinite(windows).all()
        all_windows.append(windows);tail=states[-2048:]
        returns=[float(np.max(np.linalg.norm(tail[lag:]-tail[:-lag],axis=1))) for lag in range(1,65)]
        period=next((k+1 for k,r in enumerate(returns) if r<1e-10*(1+np.linalg.norm(np.ptp(tail,axis=0)))),None)
        k=list(reference['ids']).index(ident);r=reference['tails'][k];dist=cKDTree(r).query(states[-20000:])[0]/np.linalg.norm(np.ptp(r,axis=0))
        x,m,u=.73,0.,0.;first=[];b,c,e=p
        for t in range(1,21):
            m=b*m+(1-b)*x;u=max(c*u,abs(x)+e);x-=m/((1-b**t)*u);first.append([x,m,u])
        first=np.array(first);initial_error=np.max(abs(first-source['first20'][i])/(1+abs(first)))
        row=dict(id=int(ident),parameters=p.tolist(),autonomous_step_max_scaled_residual=float(error),
            spectrum_two=windows.mean(0).tolist(),leading_windows=windows[:,0].tolist(),minimum_leading_window=float(windows[:,0].min()),
            observed_after_alignment=190000,tangent_alignment_discard=10000,reset_fraction=float(resets/190000),
            tail_period_up_to64=period,support_distances_q50_q95_q99_max=np.quantile(dist,[.5,.95,.99,1]).tolist(),
            first20_normalization_max_scaled_difference=float(initial_error))
        rows.append(row);print(ident,row['spectrum_two'],'period',period,'map residual',error,flush=True)
    np.savez_compressed(out/'windows.npz',ids=source['ids'],windows=all_windows)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows,completed=len(rows),declared=4,
        source=str(args.source.resolve()),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        reference=str(args.reference.resolve()),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        scope='The actual upstream optimizer generated every state. Analytic branch tangents evaluated along those states; no expectation of pathwise agreement between chaotic runs.'),indent=2)+'\n')
    for name in ['analyze_adamax_upstream.py','adamax_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

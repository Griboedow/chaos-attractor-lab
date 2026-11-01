"""Multiple-shooting return proposals for G=F^4; no manifold certificate."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import root
from product_return_map import iterate


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    data=dict(np.load(args.source/'trace.npz'));report=json.loads((args.source/'report.json').read_text(encoding='utf-8'))
    c=data['center'];B=data['basis'];C=np.linalg.inv(B);beta,tau=data['params'];points=data['points'];attempts=[]
    selected=sorted([r for r in report['attempts'] if r.get('root_residual',1)<1e-8 and r.get('end_distance',1)<.04],key=lambda r:r['end_distance'])[:12]
    for seed_index,seed in enumerate(selected):
        first=[c+seed['parameter']*B[:,0]]
        for k in range(seed['steps']):first.append(iterate(first[-1],beta,tau,points,1)[0])
        for tail in [8,12,16]:
            guess=list(first);end=C@(guess[-1]-c);end[0]=0
            for k in range(tail):
                end[1:]*=data['multipliers'][1:];guess.append(c+B@end)
            guess=np.array(guess);N=len(guess)-1
            def evaluate(v,jac=False):
                z=v.reshape(-1,3);f=np.zeros(3*(N+1));A=np.zeros((3*(N+1),3*(N+1)))
                for k in range(N):
                    y,J,r,flag=iterate(z[k],beta,tau,points,1)
                    assert not flag
                    f[3*k:3*k+3]=y-z[k+1];A[3*k:3*k+3,3*k:3*k+3]=J;A[3*k:3*k+3,3*k+3:3*k+6]=-np.eye(3)
                f[-3:-1]=(C@(z[0]-c))[1:];f[-1]=(C@(z[-1]-c))[0]
                A[-3:-1,:3]=C[1:];A[-1,-3:]=C[0]
                return A if jac else f
            fit=root(evaluate,guess.ravel(),jac=lambda v:evaluate(v,True),method='hybr',options={'xtol':1e-10,'maxfev':600})
            orbit=fit.x.reshape(-1,3);residual=float(np.max(abs(evaluate(fit.x))));excursion=float(np.max(np.linalg.norm(orbit-c,axis=1)));parameter=float((C@(orbit[0]-c))[0]);terminal=float(np.linalg.norm(orbit[-1]-c))
            valid=residual<1e-10 and excursion>.04 and abs(parameter)>1e-9
            row=dict(id=len(attempts),seed_index=seed_index,seed=seed,tail=tail,steps=N,solver_success=bool(fit.success),evaluations=int(fit.nfev),root_residual=residual,max_excursion=excursion,parameter=parameter,end_distance=terminal,grid_index=seed['grid_index'],accepted_numerical=valid)
            if valid:row['orbit']=orbit.tolist()
            attempts.append(row);print(json.dumps({k:v for k,v in row.items() if k not in ['seed','orbit']}),flush=True)
            result=dict(source_sha256=hashlib.sha256((args.source/'trace.npz').read_bytes()).hexdigest(),attempts=attempts,seconds=time.time()-begin,scope='Multiple shooting with linear initial unstable and terminal stable constraints. These are finite numerical window guides, not validated invariant manifolds or homoclinic orbits.')
            (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    (args.out/'trace.npz').write_bytes((args.source/'trace.npz').read_bytes())
    for name in [Path(__file__).name,'product_return_map.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

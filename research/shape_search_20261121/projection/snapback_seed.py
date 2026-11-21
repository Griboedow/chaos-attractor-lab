"""Locate regular numerical preimage chains returning to a periodic repeller."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.spatial import cKDTree
from cycle_seeds import orbit_eval
from reproduce import raw_parameters


def inverse(target,seed,raw,n):
    def fun(q):
        y,j,_=orbit_eval(q,raw,n)
        return y-target,j
    try:
        sol=root(fun,seed,jac=True,tol=1e-12)
        res,j=fun(sol.x)
        if max(abs(res))>2e-9 or not np.isfinite(res).all(): return
        return sol.x
    except (ValueError,ArithmeticError,np.linalg.LinAlgError): return


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    rec=json.loads((args.source/f'case-{args.case}.json').read_text());raw=raw_parameters(rec)
    with np.load(rec['cloud']) as f: z=f['points'].astype(float)
    certs=json.loads((args.source/'cycles'/f'cycles-{args.case}.json').read_text())['cycles']
    cert=next(c for c in certs if c['proof']['unstable_dimension']==3)
    period=cert['period']; base=np.array(cert['center'],float)
    _,_,orb=orbit_eval(base,raw,period)
    tree=cKDTree(z);rng=np.random.default_rng(371+args.case)
    for phase in range(period):
        x=orb[phase];_,mat,_=orbit_eval(x,raw,period)
        vals,vec=np.linalg.eig(mat)
        if np.max(abs(vals.imag))<1e-9:
            basis=vec.real
        else:
            real=np.flatnonzero(abs(vals.imag)<1e-9)[0];com=np.argmax(vals.imag)
            basis=np.column_stack([vec[:,real].real,vec[:,com].real,vec[:,com].imag])
        # Preserve the real and imaginary column scale of the complex eigenvector.
        inv=np.linalg.inv(basis)
        c=np.linalg.inv(mat)
        def norm(y): return np.max(abs(inv@(y-x)))
        distances,nearest=tree.query(x,k=1500)
        records=[];best=1e9
        for depth in [period,2*period,3*period]:
            picked=[]
            for j in nearest:
                if j>=depth and all(np.linalg.norm(z[j-depth]-t)>.025 for t in picked):picked.append(z[j-depth])
                if len(picked)>=45:break
            # Additional seeds are sampled in the portion of state space visited
            # by the trajectory and in the vicinity of the known repeller.
            picked.extend(z[rng.choice(len(z),80,replace=False)])
            picked.extend(x+rng.normal(0,.15,(35,3)))
            merged=[]
            for seed in picked:
                w=inverse(x,seed,raw,depth)
                if w is None or min(np.linalg.norm(w-q) for q in orb)<1e-5 or any(np.linalg.norm(w-a)<1e-6 for a in merged):continue
                merged.append(w)
                blocks=[dict(period=depth,seed=w.tolist())]
                previous=norm(w)
                for k in range(46):
                    base_seed=x+c@(w-x)
                    trial=[base_seed,x+.5*(base_seed-x),x+.2*(base_seed-x),x]
                    preds=[q for seed in trial if (q:=inverse(w,seed,raw,period)) is not None]
                    if not preds:break
                    new=min(preds,key=norm);d=norm(new)
                    if d>max(.3,2*previous):break
                    blocks.append(dict(period=period,seed=new.tolist()))
                    previous=d;w=new
                    if d<best:
                        best=d
                        if d<2e-7:
                            data=dict(case=args.case,phase=phase,period=period,base=x.tolist(),basis=basis.tolist(),distance=d,blocks=blocks)
                            (args.out/f'snap-{args.case}-{phase}.json').write_text(json.dumps(data,indent=2))
                    if d<2e-11: break
                if previous<1e-7:
                    print('found',args.case,phase,depth,len(blocks),previous,flush=True)
                    return
            records.append([depth,len(merged)])
        print('phase',phase,'best distance',best,'root counts',records,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--case',type=int,required=True)
    main(p.parse_args())

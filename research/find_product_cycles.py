"""Numerical periodic-orbit proposals for the product-alphabet map."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit
from scipy.optimize import root
from search_product_alphabet import step


@njit(cache=True)
def iterate(x,beta,tau,levels,order,period):
    J=np.eye(3);minimum=np.inf
    for k in range(period):
        x,D,r,_,_,flag=step(x,beta,tau,levels,order);minimum=min(minimum,r)
        if flag or not np.isfinite(x).all() or np.linalg.norm(x)>100:return x,J,minimum,1
        J=D@J
    return x,J,minimum,0


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time();attempts=[];cycles=[]
    with np.load(args.source) as f:beta,tau=f['params'];levels=f['levels'];order=int(f['order']);points=f['points'];z=f['states']
    for name in [Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for period in args.periods:
        distances=np.linalg.norm(z[period:]-z[:-period],axis=1);seeds=[]
        for index in np.argsort(distances):
            if all(np.linalg.norm(z[index]-z[j])>1e-3 for j in seeds):seeds.append(int(index))
            if len(seeds)==args.seeds:break
        for index in seeds:
            def function(x):
                y,J,r,flag=iterate(x,beta,tau,levels,order,period)
                return np.full(3,1e5) if flag else y-x
            def jacobian(x):return iterate(x,beta,tau,levels,order,period)[1]-np.eye(3)
            fit=root(function,z[index],jac=jacobian,method='hybr',options={'xtol':1e-11,'maxfev':250})
            y,J,r,flag=iterate(fit.x,beta,tau,levels,order,period);residual=float(np.max(abs(y-fit.x)))
            row=dict(period=period,seed_index=index,seed_return=float(distances[index]),solver_success=bool(fit.success),status=int(flag),residual=residual,evaluations=int(fit.nfev))
            if not flag and residual<1e-10:
                primitive=period
                for d in range(1,period):
                    if period%d==0 and np.linalg.norm(iterate(fit.x,beta,tau,levels,order,d)[0]-fit.x)<1e-8:primitive=d;break
                q=fit.x.copy();orbit=[]
                for k in range(primitive):orbit.append(q);q=iterate(q,beta,tau,levels,order,1)[0]
                orbit=np.array(orbit);duplicate=next((i for i,c in enumerate(cycles) if c['period']==primitive and np.min(np.linalg.norm(c['orbit']-fit.x,axis=1))<1e-7),None)
                if duplicate is None:
                    orbit=np.roll(orbit,-np.lexsort((orbit[:,2],orbit[:,1],orbit[:,0]))[0],axis=0)
                    y,J,r,flag=iterate(orbit[0],beta,tau,levels,order,primitive);values,vectors=np.linalg.eig(J);indices=np.argsort(-abs(values));values=values[indices];vectors=vectors[:,indices]
                    ident=len(cycles);cycle=dict(id=ident,period=primitive,orbit=orbit,residual=float(np.max(abs(y-orbit[0]))),minimum_radius=float(r),multipliers=[[float(v.real),float(v.imag)] for v in values],unstable_directions=int((abs(values)>1).sum()),eigenvector_condition=float(np.linalg.cond(vectors)))
                    np.savez_compressed(args.out/f'cycle-{ident}.npz',params=[beta,tau],levels=levels,order=order,points=points,orbit=orbit,monodromy=J,multipliers=values,eigenvectors=vectors);cycles.append(cycle);duplicate=ident
                row['cycle']=duplicate
            attempts.append(row)
        print(json.dumps(dict(period=period,attempts=len(attempts),cycles=len(cycles),seconds=time.time()-begin)),flush=True)
        serial=[{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in c.items()} for c in cycles]
        report=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),periods=args.periods,seeds_per_period=args.seeds,attempts=attempts,cycles=serial,seconds=time.time()-begin,scope='Floating-point periodic-root and multiplier proposals. No interval existence, invariant-set membership or entropy certificate.')
        (args.out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--periods',type=int,nargs='+',default=list(range(1,17)));p.add_argument('--seeds',type=int,default=40);main(p.parse_args())

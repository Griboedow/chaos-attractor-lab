"""Numerical periodic-point seeds for a soft-map covering search, not proofs."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit
from scipy.optimize import root
from search_soft_constellations import step


@njit(cache=True)
def iterate(x,beta,tau,points,period):
    J=np.eye(3);minimum=1e100
    for k in range(period):
        x,D,r,peak,entropy,flag=step(x,beta,points,tau);minimum=min(minimum,r)
        if flag or not np.isfinite(x).all() or np.linalg.norm(x)>100:return x,J,minimum,1
        J=D@J
    return x,J,minimum,0


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time();attempts=[];cycles=[]
    with np.load(args.source) as f:beta,tau=f['params'];points=f['points'];z=f['states']
    with np.load(args.trap) as f:origin=f['origin'];shape=f['shape'];h=float(f['step']);occupied=set(f['occupied'].tolist())
    def inside(q):
        v=np.floor((q-origin)/h).astype(np.int64)
        return bool((v>=0).all() and (v<shape).all() and int((v[0]*shape[1]+v[1])*shape[2]+v[2]) in occupied)
    for period in args.periods:
        distances=np.linalg.norm(z[period:]-z[:-period],axis=1);order=np.argsort(distances);seeds=[]
        for index in order:
            if all(np.linalg.norm(z[index]-z[j])>1e-3 for j in seeds):seeds.append(int(index))
            if len(seeds)==args.seeds:break
        for index in seeds:
            def function(x):
                y,J,axis,flag=iterate(x,beta,tau,points,period)
                if flag:return np.full(3,1e5)
                return y-x
            def jacobian(x):return iterate(x,beta,tau,points,period)[1]-np.eye(3)
            fit=root(function,z[index],jac=jacobian,method='hybr',options={'xtol':1e-11,'maxfev':250})
            y,J,axis,flag=iterate(fit.x,beta,tau,points,period);residual=float(np.max(abs(y-fit.x)))
            row=dict(period=period,seed_index=index,seed_return=float(distances[index]),solver_success=bool(fit.success),status=int(flag),residual=residual,evaluations=int(fit.nfev))
            if not flag and residual<1e-10:
                primitive=period
                for d in range(1,period):
                    if period%d==0 and np.linalg.norm(iterate(fit.x,beta,tau,points,d)[0]-fit.x)<1e-8:primitive=d;break
                q=fit.x.copy();orbit=[]
                for k in range(primitive):orbit.append(q);q=iterate(q,beta,tau,points,1)[0]
                orbit=np.array(orbit);member=all(inside(q) for q in orbit);row.update(primitive_period=primitive,in_trapping_region=member)
                duplicate=next((i for i,c in enumerate(cycles) if c['period']==primitive and np.min(np.linalg.norm(c['orbit']-fit.x,axis=1))<1e-7),None)
                if duplicate is None:
                    canonical=np.lexsort((orbit[:,2],orbit[:,1],orbit[:,0]))[0];orbit=np.roll(orbit,-canonical,axis=0)
                    y,J,axis,flag=iterate(orbit[0],beta,tau,points,primitive);values,vectors=np.linalg.eig(J);sort=np.argsort(-abs(values));values=values[sort];vectors=vectors[:,sort]
                    ident=len(cycles);cycle=dict(id=ident,period=primitive,orbit=orbit,in_trapping_region=member,residual=float(np.max(abs(y-orbit[0]))),minimum_radius=float(axis),multipliers=[[float(v.real),float(v.imag)] for v in values],unstable_directions=int((abs(values)>1).sum()),eigenvector_condition=float(np.linalg.cond(vectors)))
                    np.savez_compressed(args.out/f'cycle-{ident}.npz',params=[beta,tau],points=points,orbit=orbit,monodromy=J,multipliers=values,eigenvectors=vectors)
                    cycles.append(cycle);duplicate=ident
                row['cycle']=duplicate
            attempts.append(row)
        print(json.dumps(dict(period=period,attempts=len(attempts),cycles=len(cycles),seconds=time.time()-begin)),flush=True)
        serial=[{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in c.items()} for c in cycles]
        result=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),trap_sha256=hashlib.sha256(args.trap.read_bytes()).hexdigest(),periods=args.periods,seeds_per_period=args.seeds,attempts=attempts,cycles=serial,seconds=time.time()-begin,scope='Floating-point root and eigenvalue proposals only. No periodic-point existence, least-period, stability or entropy certificate.')
        (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'search_soft_constellations.py','search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--trap',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--periods',type=int,nargs='+',default=list(range(1,17)));p.add_argument('--seeds',type=int,default=40);main(p.parse_args())

"""Numerical periodic-orbit continuation behind the R1051 serrated-arch slice."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import root,brentq
from numba import njit
from search_raar_circles import step


@njit(cache=True)
def iterate(x,p,period):
    s=x.copy();M=np.eye(3);distance=1e100;orbit=np.empty((period,3))
    for k in range(period):
        orbit[k]=s;s,J,dist,status=step(s,p);distance=min(distance,dist)
        if status:return s,M,distance,orbit,1
        M=J@M
    return s,M,distance,orbit,0


def solve(p,period,start):
    def fun(x):
        y,M,_,_,status=iterate(x,p,period)
        return y-x,M-np.eye(3)
    solution=root(fun,start,jac=True,tol=1e-11)
    y,M,distance,orbit,status=iterate(solution.x,p,period);residual=float(np.max(np.abs(y-solution.x)))
    good=status==0 and residual<2e-10 and np.isfinite(M).all()
    divisors=[d for d in range(1,period) if period%d==0]
    distances=[float(np.linalg.norm(orbit[d]-orbit[0])) for d in divisors]
    least=min(distances,default=float('inf'))
    return dict(state=solution.x,monodromy=M,orbit=orbit,residual=residual,minaxis=distance,complete=bool(good),proper_divisors=divisors,proper_divisor_distances=distances,least_period_separation=least,multipliers=np.linalg.eigvals(M))


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    with np.load(args.source) as d:scan={k:d[k] for k in d.files}
    source=Path(__file__);(args.out/source.name).write_bytes(source.read_bytes());(args.out/'search_raar_circles.py').write_bytes(source.with_name('search_raar_circles.py').read_bytes())
    base=scan['params'][128].copy();branches=[];certificates=[]
    for period,index,startbeta,endbeta in [(7,40,base[0]-.02,base[0]),(14,80,scan['params'][70,0],scan['params'][100,0]),(28,90,scan['params'][90,0],scan['params'][105,0]),(56,95,scan['params'][95,0],scan['params'][104,0])]:
        seed=solve(scan['params'][index],period,scan['ends'][index]);assert seed['complete'] and seed['least_period_separation']>1e-5
        certificates.append(dict(label=f'stable-period-{period}',period=period,params=scan['params'][index].tolist(),state=seed['state'].tolist(),source_slice_row=index))
        x=seed['state'];rows=[]
        for beta in np.linspace(startbeta,endbeta,401):
            p=base.copy();p[0]=beta;r=solve(p,period,x)
            # Rephase the same periodic orbit to the closest point to the previous state.
            if r['complete']:
                shift=np.argmin(np.linalg.norm(r['orbit']-x,axis=1));x=r['orbit'][shift].copy();r=solve(p,period,x)
            r.update(beta=float(beta),params=p.tolist());rows.append(r)
        np.savez_compressed(args.out/f'period-{period}-branch.npz',beta=[r['beta'] for r in rows],states=[r['state'] for r in rows],multipliers=[r['multipliers'] for r in rows],residual=[r['residual'] for r in rows],minaxis=[r['minaxis'] for r in rows],complete=[r['complete'] for r in rows],least_period_separation=[r['least_period_separation'] for r in rows])
        if period==7:
            r=rows[-1]
            if r['complete'] and r['least_period_separation']>1e-5:certificates.append(dict(label='base-unstable-period-7',period=7,params=base.tolist(),state=r['state'].tolist(),source_slice_row=128))
        crossings=[]
        for left,right in zip(rows,rows[1:]):
            if not(left['complete'] and right['complete'] and min(left['least_period_separation'],right['least_period_separation'])>1e-4):continue
            dl=np.linalg.det(left['monodromy']+np.eye(3));dr=np.linalg.det(right['monodromy']+np.eye(3))
            if dl*dr<0:
                def value(beta):
                    p=base.copy();p[0]=beta;r=solve(p,period,(left['state']+right['state'])/2)
                    if not r['complete']:raise ValueError('Continuation root failed')
                    return np.linalg.det(r['monodromy']+np.eye(3))
                beta=brentq(value,left['beta'],right['beta'],xtol=5e-15);p=base.copy();p[0]=beta;r=solve(p,period,(left['state']+right['state'])/2)
                crossings.append(dict(beta=beta,multipliers=[[float(z.real),float(z.imag)] for z in r['multipliers']],residual=r['residual'],state=r['state'].tolist(),bracket=[left['beta'],right['beta']]))
        info=dict(period=period,settings=len(rows),completed=sum(r['complete'] for r in rows),least_period_preserved=sum(r['complete'] and r['least_period_separation']>1e-5 for r in rows),minus_one_crossings=crossings);branches.append(info);print(json.dumps(info),flush=True)
    report=dict(branches=branches,seconds=time.time()-begin,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),slice_source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),scope='Numerical root continuation and Floquet crossings; not a validated bifurcation, cascade accumulation or chaos proof. Every sampled branch outcome is retained.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');(args.out/'certificate-input.json').write_text(json.dumps(dict(cases=certificates),indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

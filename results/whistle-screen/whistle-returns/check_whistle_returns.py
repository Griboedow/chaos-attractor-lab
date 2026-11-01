"""Independent implicit continuation of weak and budget-limited whistle endpoints."""
import argparse,hashlib,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from whistle import normalized,dimensional_field,field_jac


def write(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def run(job):
    began=time.perf_counter();p=job['physical'];coeff,scale,omega,_=normalized(p)
    x=np.array(job['initial'])
    def section(t,q):return q[4]
    section.direction=1
    sol=solve_ivp(lambda t,q:dimensional_field(q*scale,p)/(scale*omega),(0,1000),x,
        method='Radau',jac=lambda t,q:field_jac(q,coeff)[1],rtol=2e-10,atol=2e-12,
        max_step=.05,events=section,dense_output=True)
    assert sol.success
    e=sol.t_events[0];s=sol.y_events[0];tail=s[-21:];scaler=1+np.max(abs(s),axis=0)
    error=None if len(tail)<21 else float(np.max(np.linalg.norm(np.diff(tail,axis=0)/scaler,axis=1)))
    grid=np.linspace(800,1000,10001);history=sol.sol(grid).T
    path=Path(job['output'])/f"case-{job['id']:04}.npz"
    np.savez_compressed(path,time=grid,states=history,section_times=e,section_states=s,coefficients=coeff,scale=scale,omega=omega)
    row=dict(id=job['id'],source_status=job['source_status'],events=len(e),last20_full_state_scaled_lag1_max=error,
        last20_mean_period=None if len(e)<21 else float(np.mean(np.diff(e[-21:]))),
        numerical_period1=error is not None and error<1e-7,state_extrema=np.array([history.min(0),history.max(0)]).tolist(),
        file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),wall_seconds=time.perf_counter()-began)
    write(path.with_suffix('.json'),row);return row


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source/'complete-or-budgeted.npz');meta=json.loads((args.source/'protocol.json').read_text())
    ids=sorted(set([102,134]+np.flatnonzero(source['status']!=0).tolist()))
    jobs=[dict(id=i,source_status=int(source['status'][i]),initial=source['final_states'][i].tolist(),
        physical=meta['physical_parameters'][int(source['tuple_id'][i])],output=str(out)) for i in ids]
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256((args.source/'complete-or-budgeted.npz').read_bytes()).hexdigest(),
        ids=ids,jobs=jobs,method='Independent dimensional field,Radau2e-10/2e-12,maxstep.05,1000tau from saved endpoints.',
        selection='All four step-budget exits and the two completed lanes with highest minimum window estimate.',
        return_rule='At least21 positive q=0 crossings; maximum last20 consecutive full6D return difference,scaled by1+max(abs(all section states)),<1e-7.',
        scope='A periodic continuation endpoint does not complete or retroactively classify the original GPU time interval. Finite numerical return test,not global proof.'))
    for name in ['check_whistle_returns.py','whistle.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(run,j) for j in jobs]):
            row=f.result();rows.append(row);write(out/'analysis.json',dict(completed=len(rows),declared=len(jobs),runs=sorted(rows,key=lambda x:x['id'])))
            print(row,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

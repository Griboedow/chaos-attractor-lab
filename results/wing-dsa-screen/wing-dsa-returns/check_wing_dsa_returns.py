"""Independent implicit return checks for every weak wing-screen signal."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from wing_dsa import source_equivalence,field_jac


def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def run(job):
    began=time.perf_counter();p=np.array(job['parameters']);q=np.array(job['initial'])
    f,j=field_jac(q,p);source_error=float(np.max(abs(source_equivalence(q,p)-f)/(1+abs(f))))
    assert source_error<1e-10
    rhs=lambda t,x:source_equivalence(x,p)
    warm=solve_ivp(rhs,(0,500),q,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.03)
    assert warm.success
    def section(t,q):return q[0]
    section.direction=1
    times=np.linspace(0,job['duration'],1+int(round(job['duration']/.02)))
    sol=solve_ivp(rhs,(0,job['duration']),warm.y[:,-1],method='Radau',jac=lambda t,x:field_jac(x,p)[1],
        rtol=2e-10,atol=2e-12,max_step=.03,t_eval=times,events=section)
    assert sol.success
    events=sol.y_events[0];tail=events[-100:];metrics=[]
    if len(tail)>=10:
        scale=1+np.max(abs(tail),axis=0)
        for lag in range(1,min(16,len(tail)//3)+1):
            metrics.append(dict(lag=lag,maximum_scaled_distance=float(np.max(abs(tail[lag:]-tail[:-lag])/scale))))
    close=[r['lag'] for r in metrics if r['maximum_scaled_distance']<1e-7]
    path=Path(job['output'])/f"case-{job['id']:04}.npz"
    np.savez_compressed(path,parameters=p,initial=q,time=times,states=sol.y.T,section_times=sol.t_events[0],section_states=events)
    record=dict(id=job['id'],parameters=p.tolist(),events=len(events),first_numerical_period_lag=min(close) if close else None,
        return_distances=metrics,source_error=source_error,wall_seconds=time.perf_counter()-began,
        maximum_absolute_state=float(np.max(abs(sol.y))),file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    write(path.with_suffix('.json'),record);return record


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=args.source/'complete-or-budgeted.npz';data=np.load(source);analysis=json.loads((args.source/'analysis.json').read_text())
    ids=analysis['weak_ids'];jobs=[dict(id=i,parameters=data['parameters'][i].tolist(),initial=data['states'][i,:4].tolist(),duration=args.duration,output=str(out)) for i in ids]
    alignment=[]
    for i in ids:
        f=field_jac(data['states'][i,:4],data['parameters'][i])[0];v=data['states'][i,4:]
        alignment.append(dict(id=i,one_minus_absolute_flow_cosine=float(1-abs(np.dot(v,f)/(np.linalg.norm(v)*np.linalg.norm(f))))))
    write(out/'protocol.json',dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),jobs=jobs,
        selection='All completed scout lanes with all four leading-rate windows>.001; no primary .005 candidate existed.',
        method='Saved endpoints, DOP853500tau warm, Radau800tau observation in independently transformed dimensional Eq9; exact Jacobian.',
        recurrence='Solver-root x=0 upward events, last100, all4 coordinates, maximum scaled error<1e-7, lags1..16.',
        tangent_alignment=alignment,limitation='Finite numerical returns; no global absence-of-chaos theorem.'))
    for name in ['check_wing_dsa_returns.py','wing_dsa.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(run,job) for job in jobs]):
            r=future.result();records.append(r);write(out/'analysis.json',dict(declared=len(jobs),completed=len(records),runs=records))
            print(r['id'],'events',r['events'],'lag',r['first_numerical_period_lag'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--duration',type=float,default=800.);main(p.parse_args())

"""Independent Radau return checks of declared raw-coordinate polymer signals."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from supramolecular_log import field_jac,independent_field,from_raw,to_raw


def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')


def run(job):
    began=time.perf_counter();p=np.array(job['parameters'])
    initial=np.array(job['log_initial']) if 'log_initial' in job else from_raw(np.array(job['physical_initial']),p)
    f,j=field_jac(initial,p);jj=np.column_stack([independent_field(initial+1e-25j*np.eye(3)[i],p).imag/1e-25 for i in range(3)])
    error=float(max(np.max(abs(f-independent_field(initial,p))/(1+abs(f))),np.max(abs(j-jj)/(1+abs(jj)))))
    assert error<1e-10
    rhs=lambda t,q:field_jac(q,p)[0];jac=lambda t,q:field_jac(q,p)[1]
    warm=solve_ivp(rhs,(0,job.get('warm',500)),initial,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.02)
    assert warm.success
    def section(t,q):return q[0]-1
    section.direction=1
    grid=np.linspace(0,job['duration'],1+int(round(job['duration']/.005)))
    sol=solve_ivp(rhs,(0,job['duration']),warm.y[:,-1],method='Radau',jac=jac,
        rtol=2e-10,atol=2e-12,max_step=.02,events=section,t_eval=grid)
    assert sol.success
    states=sol.y.T;events=sol.y_events[0];count=len(events);tail=events[-min(count,100):]
    metrics=[]
    if len(tail)>=10:
        scale=1+np.max(abs(tail),axis=0)
        for lag in range(1,min(16,len(tail)//3)+1):
            distances=np.max(abs(tail[lag:]-tail[:-lag])/scale,axis=1)
            metrics.append(dict(lag=lag,maximum_scaled_distance=float(distances.max()),median_scaled_distance=float(np.median(distances))))
    candidates=[r['lag'] for r in metrics if r['maximum_scaled_distance']<1e-7]
    raw=to_raw(states,p);length=2+np.exp(states[:,2]);polymer_mass=np.exp(states[:,1])*length/p[1]
    assert np.min(raw[:,:2])>=0 and np.max(raw[:,0]+polymer_mass)<=p[2]+1e-10
    path=Path(job['output'])/f"case-{job['id']:04}.npz"
    np.savez_compressed(path,time=grid,log_states=states,physical_states=raw,parameters=p,
        physical_initial=job['physical_initial'],log_initial=initial,section_times=sol.t_events[0],section_log_states=events)
    result=dict(id=job['id'],parameters=p.tolist(),events=count,first_numerical_period_lag=min(candidates) if candidates else None,
        return_distances=metrics,minimum_log_polymer_number=float(states[:,1].min()),
        minimum_polymer_number=float(np.exp(states[:,1].min())),minimum_length=float(length.min()),
        maximum_length=float(length.max()),mean_tail_return_interval=float(np.diff(sol.t_events[0][-100:]).mean()) if count>1 else None,
        source_algebra_error=error,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),file=path.name,
        wall_seconds=time.perf_counter()-began,
        limitation='Finite numerical return test, not an exact periodic-orbit proof. Approximate long-polymer closure can fail when polymer number becomes very small.')
    write(path.with_suffix('.json'),result);return result


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=args.source/('complete-or-budgeted.npz' if args.log_endpoints else 'declared-inputs.npz');d=np.load(source)
    jobs=[]
    for i in args.ids:
        if args.log_endpoints:
            index=int(np.flatnonzero(d['source_ids']==i)[0]);p=d['parameters'][index];q=d['states'][index,:3]
            jobs.append(dict(id=i,parameters=p.tolist(),physical_initial=to_raw(q,p).tolist(),log_initial=q.tolist(),warm=50.,duration=args.duration,output=str(out)))
        else:
            jobs.append(dict(id=i,parameters=d['parameters'][i].tolist(),physical_initial=d['initial'][i].tolist(),duration=args.duration,output=str(out)))
    write(out/'protocol.json',dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),jobs=jobs,
        selection=('All four step-budget exits from the long log-coordinate continuation.' if args.log_endpoints else 'Explicit source control and strongest raw-coordinate signals, including cases near the physical approximation limit; selection before log-coordinate recheck outcome.'),
        method=('Saved log endpoints and DOP853 warm50tau; ' if args.log_endpoints else 'Fresh physical starts and DOP853 warm500tau; ')+'warm at2e-11/.02; independent Radau observation at2e-10/.02 with exact log-coordinate Jacobian.',
        section='x=1 upward, solver root events, all returned coordinates retained.',
        recurrence='Last100 section points; max coordinate difference scaled by1+max absolute coordinate <1e-7, at least10 events; lags1..16.',
        caveat='No asymptotic chaos or novelty decision from finite-time signals.'))
    for name in ['check_supramolecular_returns.py','supramolecular_log.py','supramolecular.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            r=future.result();records.append(r);write(out/'analysis.json',dict(declared=len(jobs),completed=len(records),runs=records))
            print(r['id'],'events',r['events'],'lag',r['first_numerical_period_lag'],'minlogn',r['minimum_log_polymer_number'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',default=[0,1958,3170,3316,415]);p.add_argument('--duration',type=float,default=200.)
    p.add_argument('--log-endpoints',action='store_true')
    main(p.parse_args())

"""Independent Radau state-only returns for explicitly selected chemostat lanes."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from crossfeeding import field_jac,balance_weights,balance_in
from search_crossfeeding import write


def run(job):
    began=time.perf_counter();p=np.array(job['parameters']);x=np.array(job['initial'])
    def event(t,y):return y[0]-job['section_level']
    event.direction=1
    grid=np.linspace(0,job['duration'],20001)
    sol=solve_ivp(lambda t,y:field_jac(y,p)[0],(0,job['duration']),x,
        jac=lambda t,y:field_jac(y,p)[1],method='Radau',rtol=2e-10,atol=2e-12,max_step=.05,
        t_eval=grid,events=event)
    assert sol.success
    events=np.asarray(sol.y_events[0]).reshape(-1,5);tail=events[-max(100,len(events)//2):]
    scale=1+float(np.max(abs(tail))) if len(tail) else 1.;tests=[]
    for lag in range(1,min(33,len(tail)//4)):
        distance=float(np.max(np.linalg.norm(tail[lag:]-tail[:-lag],axis=1)))
        tests.append(dict(lag=lag,maximum_distance=distance,scaled_distance=distance/scale))
    repeated=[t for t in tests if t['scaled_distance']<1e-7]
    weights=balance_weights(p);supply=balance_in(p)
    budget_error=float(max(abs(weights@sol.y-(supply+(weights@x-supply)*np.exp(-grid)))))
    assert sol.y.min()>-1e-8 and budget_error<1e-6
    path=Path(job['output'])/f"case-{job['id']:05}.npz"
    np.savez_compressed(path,parameters=p,time_tau=grid,states=sol.y.T,section_times=sol.t_events[0],section_states=events)
    result=dict(id=job['id'],events=len(events),return_tests=tests,first_numerical_period_lag=repeated[0]['lag'] if repeated else None,
        budget_max_error=budget_error,minimum=sol.y.min(axis=1).tolist(),maximum=sol.y.max(axis=1).tolist(),
        file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),wall_seconds=time.perf_counter()-began,
        interpretation='Finite numerical recurrence screen; neither a periodic-orbit proof nor an absence-of-chaos theorem.')
    write(path.with_suffix('.json'),result);return result


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);data=np.load(args.source)
    jobs=[]
    for i in args.ids:
        level=float(data['stats'][i,:2].mean())
        assert data['status'][i]==0 and np.isfinite(level)
        jobs.append(dict(id=i,parameters=data['parameters'][i].tolist(),initial=data['states'][i,:5].tolist(),
                         section_level=level,duration=args.duration,output=str(out)))
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),jobs=jobs,
        selection='Three largest mean leading rates among the complete weak-positive pilot lanes, plus the source control.',
        method='Independent Radau state integration with analytic piecewise Jacobian; rtol2e-10,atol2e-12,maxstep.05tau.',
        section='Upward N1 crossing of its pilot observation midpoint; all five state components retained.',
        recurrence='Lags1..32 on the last max(100,half) events; scaled maximum distance<1e-7; at least four lag blocks.',
        time='tau=D*t; author t in hours',scope='Failure to recur is inconclusive. No Lyapunov estimate is inferred from this screen.'))
    for name in ['crossfeeding.py','screen_crossfeeding_returns.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            r=future.result();records.append(r);write(out/'analysis.json',dict(declared=len(jobs),completed=len(records),runs=sorted(records,key=lambda x:x['id'])))
            print(r['id'],'events',r['events'],'return_lag',r['first_numerical_period_lag'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',default=[0,2327,4460,6143]);p.add_argument('--duration',type=float,default=400.)
    main(p.parse_args())

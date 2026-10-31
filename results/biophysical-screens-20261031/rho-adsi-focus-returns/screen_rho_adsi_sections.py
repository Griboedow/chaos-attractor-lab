"""Independent state-only return screen; no Lyapunov estimates are inferred here."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from rho_adsi_kinetics import parameterize,equilibria,independent_field
from search_rho_adsi import write


@njit(cache=True)
def rhs(t,x,k):
    active,inactive,actin=x
    a=k[0]+k[1]*active**3/(1+k[2]*active*active)
    b=k[3]+k[4]*actin;rx=a*inactive-b*active
    return np.array([rx,-rx+k[5]-k[6]*inactive,
        k[7]+k[8]*active*active/(1+k[9]*active*active)-k[10]*actin])/.025


@njit(cache=True)
def jac(t,x,k):
    active,inactive,actin=x;den=1+k[2]*active*active
    a=k[0]+k[1]*active**3/den;da=k[1]*active*active*(3+k[2]*active*active)/den**2
    b=k[3]+k[4]*actin
    return np.array([[da*inactive-b,a,-k[4]*active],[-da*inactive+b,-a-k[6],k[4]*active],
        [2*k[8]*active/(1+k[9]*active*active)**2,0.,-k[10]]])/.025


def run(job):
    began=time.perf_counter();k=np.array(job['kinetics']);x=np.array(job['initial'])
    roots=equilibria(k);level=float(np.median([r['state'][0] for r in roots]))
    def section(t,y,k):return y[0]-level
    section.direction=1
    test=independent_field(x,k);jtest=np.column_stack([independent_field(x+1e-25j*np.eye(3)[i],k).imag/1e-25 for i in range(3)])
    algebra=float(max(np.max(abs(rhs(0,x,k)-test)/(1+abs(test))),np.max(abs(jac(0,x,k)-jtest)/(1+abs(jtest)))))
    assert algebra<1e-11
    sol=solve_ivp(rhs,(0,job['duration']),x,args=(k,),jac=jac,method='Radau',rtol=2e-10,atol=2e-12,
                  max_step=.1,events=section,dense_output=True)
    assert sol.success
    event_times=sol.t_events[0];event_states=np.asarray(sol.y_events[0]).reshape(-1,3)
    grid=np.linspace(0,job['duration'],20001);states=sol.sol(grid).T
    tail=event_states[-max(100,len(event_states)//2):];scale=1+float(np.max(abs(tail))) if len(tail) else 1.
    recurrence=[]
    for lag in range(1,min(33,len(tail)//4)):
        error=float(np.max(np.linalg.norm(tail[lag:]-tail[:-lag],axis=1)))
        recurrence.append(dict(lag=lag,maximum_distance=error,scaled_distance=error/scale))
    periodic=[r for r in recurrence if r['scaled_distance']<1e-7]
    out=Path(job['output']);path=out/f"case-{job['id']:05}.npz"
    np.savez_compressed(path,kinetics=k,time_tau=grid,states=states,section_level=level,
                        section_times=event_times,section_states=event_states)
    record=dict(id=job['id'],events=len(event_states),section_level=level,return_tests=recurrence,
        first_numerical_period_lag=periodic[0]['lag'] if periodic else None,algebra_relative_error=algebra,
        minimum=states.min(axis=0).tolist(),maximum=states.max(axis=0).tolist(),
        wall_seconds=time.perf_counter()-began,file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        interpretation='Finite numerical return screen, not a periodic-orbit proof or an absence-of-chaos theorem.')
    write(path.with_suffix('.json'),record);return record


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);data=np.load(args.source)
    ids=args.ids or np.flatnonzero((data['status']==0)&(data['windows'].min(axis=1)>.001)).tolist()
    jobs=[]
    for i in ids:
        assert data['status'][i]==0
        p=data['parameters'][i]
        k=parameterize(np.array([[p[0],p[1]/p[0],p[2],p[3],1,1,1,1]]))[0] if len(p)==4 else p
        jobs.append(dict(id=i,kinetics=k.tolist(),initial=data['states'][i,:3].tolist(),duration=args.duration,output=str(out)))
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        jobs=jobs,method='Independent Radau state integration with analytic Jacobian;rtol2e-10,atol2e-12,maxstep.1tau.',
        selection='All complete exploratory weak lanes with all4 windows>.001, unless explicit IDs supplied.',
        time='tau=.025*t, physical t in seconds',diagnostic='UpwardRT section at median positive equilibrium RT; recurrence lags1..32 across last max(100,half) events; scaled tolerance1e-7.',
        scope='Failure to recur is inconclusive. No positive-exponent assertion from this screen.'))
    for name in ['screen_rho_adsi_sections.py','rho_adsi_kinetics.py','rho_adsi.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(run,job) for job in jobs]):
            r=future.result();records.append(r)
            write(out/'analysis.json',dict(declared=len(jobs),completed=len(records),runs=sorted(records,key=lambda r:r['id'])))
            print(r['id'],'events',r['events'],'return_lag',r['first_numerical_period_lag'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+')
    p.add_argument('--duration',type=float,default=1000.)
    main(p.parse_args())

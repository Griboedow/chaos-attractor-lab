"""Declared local N=3 geometry survey with full finite-time spectra.

One carried historical start at each parameter point is not a basin survey.
All outcomes are retained, including regular and weakly resolved cases.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from motor_rows import field_jac, rhs

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,record):
    path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def augmented(z,p):
    rho,r,b=p
    c=z[:3];s=z[3:6];x=np.array([z[6],z[7],-z[6]-z[7]])
    v=-rho*x+r*(s-s.mean())
    f=np.concatenate((1-c-v*s,b-s+v*c,v[:2]))
    m=z[8:72].reshape(8,8)
    dx=np.vstack((m[6:8],-m[6:8].sum(axis=0)))
    dv=-rho*dx+r*(m[3:6]-m[3:6].mean(axis=0))
    dm=np.vstack((-m[:3]-v[:,None]*m[3:6]-s[:,None]*dv,
                   -m[3:6]+v[:,None]*m[:3]+c[:,None]*dv,dv[:2]))
    return np.concatenate((f,dm.ravel(),[-6-2*rho+(2*r/3)*c.sum()]))


def run(job):
    began=time.perf_counter();p=np.array(job['parameters']);x=np.array(job['initial'])
    out=Path(job['output']);basis=np.eye(8);logs=np.zeros(8);trace=0.;elapsed=0.
    test=np.r_[x,basis.ravel(),0.];f,j=field_jac(x,p,3)
    got=augmented(test,p)
    algebra=max(float(np.max(abs(got[:8]-f))),float(np.max(abs(got[8:72].reshape(8,8)-j))),abs(got[-1]-np.trace(j)))
    assert algebra<1e-12
    times=[];states=[];section_times=[];section_states=[];windows=[]
    previous=np.zeros(8);previous_trace=0.
    warm=job['warm'];duration=job['duration'];interval=.5
    def section(t,z):return z[6]
    section.direction=1
    while elapsed<warm+duration-1e-9:
        end=min(elapsed+interval,warm+duration)
        sample=np.linspace(elapsed,end,26)
        sol=solve_ivp(lambda t,z:augmented(z,p),(elapsed,end),np.r_[x,basis.ravel(),0.],
                      method='DOP853',rtol=2e-10,atol=2e-12,max_step=.04,t_eval=sample,events=section)
        if not sol.success:raise RuntimeError(sol.message)
        x=sol.y[:8,-1];basis,tri=np.linalg.qr(sol.y[8:72,-1].reshape(8,8))
        increments=np.log(abs(np.diag(tri)))
        if elapsed>=warm-1e-9:
            logs+=increments;trace+=sol.y[-1,-1]
            times.extend(sample[1:]);states.extend(sol.y[:8,1:].T)
            section_times.extend(sol.t_events[0])
            section_states.extend(np.asarray(sol.y_events[0]).reshape(-1,73)[:,:8])
        elapsed=end
        if elapsed>warm and abs((elapsed-warm)%250)<1e-8:
            windows.append(dict(end=elapsed,spectrum=((logs-previous)/250).tolist(),
                                mean_divergence=(trace-previous_trace)/250))
            previous=logs.copy();previous_trace=trace
    xs=np.array(states);c=xs[:,:3];s=xs[:,3:6];X=np.column_stack((xs[:,6:],-xs[:,6:].sum(axis=1)))
    distances=np.column_stack([np.sqrt((c[:,a]-c[:,b])**2+(s[:,a]-s[:,b])**2+(X[:,a]-X[:,b])**2)
                               for a,b in [(0,1),(0,2),(1,2)]])
    path=out/f"case-{job['id']:02}.npz"
    np.savez_compressed(path,time=times,states=xs,parameters=p,section_time=section_times,
                        section_states=np.array(section_states).reshape(-1,8),final_state=x,
                        final_basis=basis,log_sums=logs,integrated_trace=trace)
    result=dict(id=job['id'],parameters=p.tolist(),file=path.name,sha256=sha(path),
                spectrum=(logs/duration).tolist(),windows=windows,algebra_error=algebra,
                trace_residual_per_tau=float(abs(logs.sum()-trace)/duration),
                max_disk_radius_squared=float(np.max(c*c+s*s)),
                pair_rms=np.sqrt(np.mean(distances**2,axis=0)).tolist(),
                sections=len(section_times),samples=len(times),wall_seconds=time.perf_counter()-began,
                status='finite-time numerical observation; no attraction or novelty proof')
    assert np.isfinite(xs).all() and result['max_disk_radius_squared']<=2+1e-8
    assert result['trace_residual_per_tau']<1e-6
    write(out/f"case-{job['id']:02}.json",result)
    return result


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    src=ROOT/'results/motor-rows/cpu-candidate-N3-783.npz';z=np.load(src)
    rho0,r0,b=z['parameters'];q0=r0/(rho0+1)
    jobs=[]
    for rho in [3.7,float(rho0),5.6]:
        for q in [1.13,float(q0),1.5]:
            jobs.append(dict(id=len(jobs),parameters=[rho,q*(rho+1),1.],initial=z['final_state'].tolist(),
                             output=str(out),warm=200.,duration=1000.))
    protocol=dict(source=src.relative_to(ROOT).as_posix(),source_sha256=sha(src),script_sha256=sha(__file__),
                  source_field_sha256=sha(ROOT/'research/motor_rows.py'),jobs=jobs,
                  selection='Cartesian 3 by 3 grid declared before integration; one common historical endpoint.',
                  integrator='DOP853 full 8D QR; QR 0.5; maxstep .04; rtol2e-10,atol2e-12',
                  sampling='Uniform .02 original tau throughout the observed interval.',
                  scope='Dimensionless mathematical family at b=1, not the one-parameter physical activity path.')
    write(out/'protocol.json',protocol);(out/'source.py').write_bytes(Path(__file__).read_bytes())
    (out/'motor_rows.py').write_bytes((ROOT/'research/motor_rows.py').read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run,j) for j in jobs]
        for future in as_completed(futures):
            item=future.result();records.append(item)
            write(out/'analysis.json',dict(completed=len(records),declared=len(jobs),runs=sorted(records,key=lambda x:x['id'])))
            print('case',item['id'],'spectrum',np.round(item['spectrum'][:3],6),'sections',item['sections'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--workers',type=int,default=3);main(p.parse_args())

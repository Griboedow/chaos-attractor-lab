"""Full-spectrum physical-geometry follow-up of stratified motor scouts."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from motor_rows import field_jac


@njit(cache=True)
def augmented(z,p):
    rho,r,b=p
    answer=np.empty(73)
    mean_s=(z[3]+z[4]+z[5])/3
    v=np.empty(3)
    for i in range(3):
        x=z[6+i] if i<2 else -z[6]-z[7]
        v[i]=-rho*x+r*(z[3+i]-mean_s)
        answer[i]=1-z[i]-v[i]*z[3+i]
        answer[3+i]=b-z[3+i]+v[i]*z[i]
        if i<2:answer[6+i]=v[i]
    for col in range(8):
        mean_ds=(z[8+3*8+col]+z[8+4*8+col]+z[8+5*8+col])/3
        for i in range(3):
            dc=z[8+i*8+col]
            ds=z[8+(3+i)*8+col]
            dx=z[8+(6+i)*8+col] if i<2 else -z[8+6*8+col]-z[8+7*8+col]
            dv=-rho*dx+r*(ds-mean_ds)
            answer[8+i*8+col]=-dc-v[i]*ds-z[3+i]*dv
            answer[8+(3+i)*8+col]=-ds+v[i]*dc+z[i]*dv
            if i<2:answer[8+(6+i)*8+col]=dv
    answer[-1]=-6-2*rho+(2*r/3)*(z[0]+z[1]+z[2])
    return answer


def write(p,x):
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')


def run(job):
    began=time.perf_counter()
    p=np.array(job['parameters']);x=np.array(job['initial']);basis=np.eye(8)
    reference,j=field_jac(x,p,3)
    got=augmented(np.r_[x,basis.ravel(),0.],p)
    error=max(float(max(abs(got[:8]-reference))),float(np.max(abs(got[8:72].reshape(8,8)-j))))
    assert error<1e-11
    logs=np.zeros(8);trace=0.;states=[];times=[];windows=[]
    section_times=[];section_states=[];previous=np.zeros(8);previous_trace=0.
    warm_steps=int(round(job['warm']/.5));observed_steps=int(round(job['duration']/.5))
    assert observed_steps%4==0
    def section(t,z):return z[6]
    section.direction=1
    for it in range(warm_steps+observed_steps):
        start,end=it*.5,(it+1)*.5;grid=np.linspace(start,end,26)
        sol=solve_ivp(lambda t,z:augmented(z,p),(start,end),np.r_[x,basis.ravel(),0.],
                      method='DOP853',rtol=2e-10,atol=2e-12,max_step=.02,t_eval=grid,events=section)
        if not sol.success:raise RuntimeError(sol.message)
        x=sol.y[:8,-1];basis,tri=np.linalg.qr(sol.y[8:72,-1].reshape(8,8))
        if it>=warm_steps:
            logs+=np.log(abs(np.diag(tri)));trace+=sol.y[-1,-1]
            times.extend(grid[1:]);states.extend(sol.y[:8,1:].T)
            section_times.extend(sol.t_events[0]);section_states.extend(np.asarray(sol.y_events[0]).reshape(-1,73)[:,:8])
        if it+1>warm_steps and (it+1-warm_steps)%(observed_steps//4)==0:
            windows.append(dict(end=end,spectrum=((logs-previous)/(job['duration']/4)).tolist()))
            previous=logs.copy();previous_trace=trace
    states=np.array(states)
    path=Path(job['output'])/f"case-{job['id']:04}.npz"
    np.savez_compressed(path,time=times,states=states,parameters=p,final_state=x,final_basis=basis,
                        log_sums=logs,integrated_trace=trace,section_times=section_times,
                        section_states=np.asarray(section_states).reshape(-1,8))
    spectrum=logs/job['duration'];radius=float(np.max(states[:,:3]**2+states[:,3:6]**2))
    trace_error=float(abs(logs.sum()-trace)/job['duration'])
    assert radius<=1+p[2]**2+1e-7 and trace_error<1e-6
    result=dict(id=job['id'],parameters=p.tolist(),spectrum=spectrum.tolist(),windows=windows,
                algebra_error=error,trace_residual_per_tau=trace_error,maximum_disk_radius_squared=radius,
                samples=len(states),sections=len(section_times),file=path.name,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),wall_seconds=time.perf_counter()-began)
    write(path.with_suffix('.json'),result)
    return result


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);p=source['parameters'];good=(source['status']==0)&(source['lle_windows'].min(axis=1)>.01)
    selected=[]
    # For each asymmetry, choose three rate quantiles among unique retained tuples.
    # This spans weak and strong signals instead of selecting only a largest exponent.
    for b in np.unique(p[:,2]):
        ids=np.flatnonzero(good&(p[:,2]==b))
        tuples={}
        for i in ids:
            tuples.setdefault(int(source['tuple_id'][i]),int(i))
        ids=np.array(list(tuples.values()),dtype=int)
        if not len(ids):continue
        ids=ids[np.argsort(source['lle'][ids])]
        selected.extend(ids[np.unique(np.rint(np.linspace(0,len(ids)-1,3)).astype(int))].tolist())
    jobs=[dict(id=int(i),parameters=p[i].tolist(),initial=source['final_states'][i].tolist(),
               warm=200.,duration=args.duration,output=str(out)) for i in sorted(selected)]
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        selection='For every b: first retained start per parameter tuple, sorted by leading rate, select endpoints and median. All four scout windows >.01.',
        selected=sorted(selected),jobs=jobs,method='DOP853 full 8D QR .5; rtol2e-10,atol2e-12,maxstep.02; uniform .02 samples',
        warning='Finite-time spectra and visual screening; no novel attractor claim.'))
    for name in ['check_motor_asymmetry.py','motor_rows.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            record=future.result();records.append(record)
            write(out/'analysis.json',dict(completed=len(records),declared=len(jobs),runs=sorted(records,key=lambda a:a['id'])))
            print(record['id'],np.round(record['spectrum'][:3],6),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=3)
    p.add_argument('--duration',type=float,default=1200.)
    main(p.parse_args())

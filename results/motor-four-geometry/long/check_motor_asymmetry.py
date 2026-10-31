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
def augmented(z,p,n=3):
    rho,r,b=p
    d=3*n-1;answer=np.empty(d+d*d+1)
    mean_s=np.sum(z[n:2*n])/n
    v=np.empty(n)
    for i in range(n):
        x=z[2*n+i] if i<n-1 else -np.sum(z[2*n:d])
        v[i]=-rho*x+r*(z[n+i]-mean_s)
        answer[i]=1-z[i]-v[i]*z[n+i]
        answer[n+i]=b-z[n+i]+v[i]*z[i]
        if i<n-1:answer[2*n+i]=v[i]
    for col in range(d):
        mean_ds=0.
        for j in range(n):mean_ds+=z[d+(n+j)*d+col]/n
        for i in range(n):
            dc=z[d+i*d+col]
            ds=z[d+(n+i)*d+col]
            dx=0.
            if i<n-1:dx=z[d+(2*n+i)*d+col]
            else:
                for j in range(n-1):dx-=z[d+(2*n+j)*d+col]
            dv=-rho*dx+r*(ds-mean_ds)
            answer[d+i*d+col]=-dc-v[i]*ds-z[n+i]*dv
            answer[d+(n+i)*d+col]=-ds+v[i]*dc+z[i]*dv
            if i<n-1:answer[d+(2*n+i)*d+col]=dv
    answer[-1]=-2*n-(n-1)*rho+r*(1-1/n)*np.sum(z[:n])
    return answer


def write(p,x):
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')


def run(job):
    began=time.perf_counter()
    p=np.array(job['parameters']);x=np.array(job['initial']);d=len(x);n=(d+1)//3;basis=np.eye(d)
    assert d==3*n-1
    reference,j=field_jac(x,p,n)
    got=augmented(np.r_[x,basis.ravel(),0.],p,n)
    error=max(float(max(abs(got[:d]-reference))),float(np.max(abs(got[d:-1].reshape(d,d)-j))))
    assert error<1e-11
    logs=np.zeros(d);trace=0.;states=[];times=[];windows=[]
    section_times=[];section_states=[];previous=np.zeros(d);previous_trace=0.
    warm_steps=int(round(job['warm']/.5));observed_steps=int(round(job['duration']/.5))
    assert observed_steps%4==0
    def section(t,z):return z[2*n]
    section.direction=1
    for it in range(warm_steps+observed_steps):
        start,end=it*.5,(it+1)*.5;grid=np.linspace(start,end,26)
        sol=solve_ivp(lambda t,z:augmented(z,p,n),(start,end),np.r_[x,basis.ravel(),0.],
                      method='DOP853',rtol=job['rtol'],atol=job['rtol']*.01,max_step=job['max_step'],t_eval=grid,events=section)
        if not sol.success:raise RuntimeError(sol.message)
        x=sol.y[:d,-1];basis,tri=np.linalg.qr(sol.y[d:-1,-1].reshape(d,d))
        if it>=warm_steps:
            logs+=np.log(abs(np.diag(tri)));trace+=sol.y[-1,-1]
            times.extend(grid[1:]);states.extend(sol.y[:d,1:].T)
            section_times.extend(sol.t_events[0]);section_states.extend(np.asarray(sol.y_events[0]).reshape(-1,d+d*d+1)[:,:d])
        if it+1>warm_steps and (it+1-warm_steps)%(observed_steps//4)==0:
            windows.append(dict(end=end,spectrum=((logs-previous)/(job['duration']/4)).tolist()))
            previous=logs.copy();previous_trace=trace
    states=np.array(states)
    path=Path(job['output'])/f"case-{job['id']:04}.npz"
    np.savez_compressed(path,time=times,states=states,parameters=p,final_state=x,final_basis=basis,
                        log_sums=logs,integrated_trace=trace,section_times=section_times,
                        section_states=np.asarray(section_states).reshape(-1,d))
    spectrum=logs/job['duration'];radius=float(np.max(states[:,:n]**2+states[:,n:2*n]**2))
    trace_error=float(abs(logs.sum()-trace)/job['duration'])
    assert radius<=1+p[2]**2+1e-7 and trace_error<1e-6
    result=dict(id=job['id'],N=n,dimension=d,parameters=p.tolist(),spectrum=spectrum.tolist(),windows=windows,
                algebra_error=error,trace_residual_per_tau=trace_error,maximum_disk_radius_squared=radius,
                samples=len(states),sections=len(section_times),file=path.name,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),wall_seconds=time.perf_counter()-began)
    write(path.with_suffix('.json'),result)
    return result


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);p=source['parameters'];good=(source['status']==0)&(source['lle_windows'].min(axis=1)>.01)
    d=source['final_states'].shape[1];n=(d+1)//3
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
    selection='For every b: first retained start per parameter tuple, sorted by leading rate, select endpoints and median. All four scout windows >.01.'
    if args.ids is not None:
        selected=args.ids
        assert all(0<=i<len(p) and source['status'][i]==0 for i in selected)
        selection='Explicit follow-up IDs, selected after inspection of the first full-spectrum geometry gallery.'
    jobs=[dict(id=int(i),parameters=p[i].tolist(),initial=source['final_states'][i].tolist(),
               warm=200.,duration=args.duration,output=str(out),rtol=args.rtol,max_step=args.max_step) for i in sorted(selected)]
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        selection=selection,
        selected=sorted(selected),jobs=jobs,N=n,dimension=d,method=f'DOP853 full {d}D QR .5; rtol{args.rtol},atol{args.rtol*.01},maxstep{args.max_step}; uniform .02 samples',
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
    p.add_argument('--ids',type=int,nargs='+')
    p.add_argument('--rtol',type=float,default=2e-10)
    p.add_argument('--max-step',type=float,default=.02)
    main(p.parse_args())

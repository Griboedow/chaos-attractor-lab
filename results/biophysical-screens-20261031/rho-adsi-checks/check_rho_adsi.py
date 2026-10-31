"""Independent full-spectrum and physical-coordinate follow-up of ADSI scouts."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from rho_adsi import field_jac,independent_field,equilibria


@njit(cache=True)
def augmented(z,p):
    x,y,f=z[:3];gef,gap,speed,exchange=p
    a=.00625+gef*.3125*x*x*x/(1+x*x)
    da=gef*.3125*x*x*(3+x*x)/(1+x*x)**2
    b=.0625+.05625*gap*f;rx=a*y-b*x
    j=np.array([[da*y-b,a,-.05625*gap*x],[-da*y+b,-a-exchange*.0208333,.05625*gap*x],
                [speed*.2812*x/(1+.25*x*x)**2,0.,-speed*.025]])/.025
    answer=np.empty(13);answer[0]=rx/.025;answer[1]=(-rx+exchange*(.0625-.0208333*y))/.025
    answer[2]=speed*(.001875+.1406*x*x/(1+.25*x*x)-.025*f)/.025
    answer[3:12]=(j@z[3:12].reshape((3,3))).ravel();answer[12]=np.trace(j)
    return answer


def write(path,x):
    path.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')


def run(job):
    began=time.perf_counter();out=Path(job['output']);p=np.array(job['parameters']);x=np.array(job['initial']);q=np.eye(3)
    got=augmented(np.r_[x,q.ravel(),0.],p)
    jf=np.column_stack([independent_field(x+1e-25j*np.eye(3)[i],p).imag/1e-25 for i in range(3)])
    reference=np.r_[independent_field(x,p),jf.ravel(),np.trace(jf)]
    algebra=float(np.max(abs(got-reference)/(1+abs(reference))));assert algebra<1e-11
    short={}
    for method in ['DOP853','Radau']:
        sol=solve_ivp(lambda t,y:independent_field(y,p),(0,2),x,method=method,rtol=2e-11,atol=2e-13,max_step=.005)
        assert sol.success;short[method]=sol.y[:,-1]
    discrepancy=float(np.max(abs(short['DOP853']-short['Radau'])/(1+abs(short['Radau']))));assert discrepancy<1e-7
    logs=np.zeros(3);trace=0.;previous=np.zeros(3);states=[];times=[];windows=[]
    section_states=[];section_times=[];section_level=float(np.median([e['state'][0] for e in equilibria(p)]))
    def event(t,z):return z[0]-section_level
    event.direction=1
    for i in range(int(job['warm']+job['duration'])):
        grid=np.linspace(i,i+1,21)
        sol=solve_ivp(lambda t,z:augmented(z,p),(i,i+1),np.r_[x,q.ravel(),0.],method='DOP853',
                      rtol=2e-10,atol=2e-12,max_step=.025,t_eval=grid,events=event)
        assert sol.success
        x=sol.y[:3,-1];q,r=np.linalg.qr(sol.y[3:12,-1].reshape(3,3))
        if i>=job['warm']:
            logs+=np.log(abs(np.diag(r)));trace+=sol.y[-1,-1]
            states.extend(sol.y[:3,1:].T);times.extend(grid[1:])
            section_times.extend(sol.t_events[0]);section_states.extend(np.asarray(sol.y_events[0]).reshape(-1,13)[:,:3])
            if (i+1-job['warm'])%(job['duration']/4)==0:
                windows.append(dict(end=i+1,spectrum=((logs-previous)/(job['duration']/4)).tolist()));previous=logs.copy()
    states=np.array(states);residual=float(abs(logs.sum()-trace)/job['duration']);assert states.min()>-1e-8 and residual<1e-6
    path=out/f"case-{job['id']:04}.npz"
    np.savez_compressed(path,time=times,states=states,parameters=p,log_sums=logs,integrated_trace=trace,
        final_state=x,final_basis=q,section_level=section_level,section_times=section_times,section_states=np.asarray(section_states).reshape(-1,3))
    record=dict(id=job['id'],parameters=p.tolist(),spectrum=(logs/job['duration']).tolist(),windows=windows,
        minimum=states.min(axis=0).tolist(),maximum=states.max(axis=0).tolist(),samples=len(states),sections=len(section_times),
        section_level=section_level,algebra_relative_error=algebra,short_solver_relative_error=discrepancy,
        trace_residual_per_tau=residual,wall_seconds=time.perf_counter()-began,file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    write(path.with_suffix('.json'),record);return record


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);source=np.load(args.source)
    assert args.duration%4==0
    jobs=[dict(id=i,parameters=source['parameters'][i].tolist(),initial=source['states'][i,:3].tolist(),
               warm=500.,duration=args.duration,output=str(out)) for i in args.ids]
    assert all(source['status'][i]==0 for i in args.ids)
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        jobs=jobs,selection='Exploratory follow-up of small, comparatively persistent positive finite-time rates, below the preregistered pilot retention threshold.',
        method='Full3D QR every1tau; DOP853 rtol2e-10,atol2e-12,maxstep.025; samples .05tau',
        time='tau=.025*t, with physical t in seconds; rates below are per tau',
        warning='Finite-time diagnostic; no attraction, uniqueness or priority proof.'))
    for name in ['check_rho_adsi.py','rho_adsi.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(run,job) for job in jobs]):
            result=future.result();records.append(result)
            write(out/'analysis.json',dict(declared=len(jobs),completed=len(records),runs=sorted(records,key=lambda r:r['id'])))
            print(result['id'],result['spectrum'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',default=[429,1612]);p.add_argument('--duration',type=int,default=8000)
    main(p.parse_args())

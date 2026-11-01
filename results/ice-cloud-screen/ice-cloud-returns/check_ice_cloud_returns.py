"""Independent returns for both starts at every weak ice-cloud tuple."""
import argparse,hashlib,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from ice_cloud import coefficients,field_jac,physical_field,to_physical


def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')


def run(job):
    began=time.perf_counter();p=np.array(job['parameters']);scale=np.array(job['scale']);omega=job['omega'];physical=job['physical']
    c=coefficients(**physical);initial=np.array(job['initial'])
    def fun(t,z):
        x=to_physical(z,scale);return physical_field(x,c)/(np.r_[x[:2],1.]*omega)
    def event(t,z):return z[2]-p[-1]
    event.direction=1
    sol=solve_ivp(fun,(0,400),initial,method='Radau',jac=lambda t,z:field_jac(z,p)[1],
        rtol=2e-11,atol=2e-13,max_step=.03,events=event,dense_output=True)
    assert sol.success
    times=sol.t_events[0];sections=sol.y_events[0];tail=sections[-51:]
    err=None if len(tail)<51 else float(np.max(np.linalg.norm(np.diff(tail,axis=0),axis=1)))
    grid=np.linspace(300,400,10001);z=sol.sol(grid).T;x=to_physical(z,scale)
    density=physical['pressure']/(287.058*physical['T'])
    radius=(3*x[:,1]/(4*np.pi*810*x[:,0]))**(1/3)
    path=Path(job['output'])/f"case-{job['id']:04}.npz"
    np.savez_compressed(path,time=grid,log_states=z,physical_states=x,section_times=times,section_states=sections,
        normalized_parameters=p,physical_coefficients=c,scale=scale,omega=omega)
    record=dict(id=job['id'],events=len(times),last50_full_log_state_lag1_max=err,numerical_period1=err is not None and err<1e-7,
        period_xi=None if len(tail)<51 else float(np.mean(np.diff(times[-51:]))),
        period_seconds=None if len(tail)<51 else float(np.mean(np.diff(times[-51:]))/omega),
        sampled_max_crystal_number_to_fixed_aerosol_number=float(np.max(x[:,0]*density/physical['aerosols'])),
        sampled_mean_mass_radius_m_range=[float(radius.min()),float(radius.max())],
        log_state_extrema=np.array([z.min(0),z.max(0)]).tolist(),file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        wall_seconds=time.perf_counter()-began)
    write(path.with_suffix('.json'),record);return record


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    d=np.load(args.source/'complete-or-budgeted.npz');inputs=np.load(args.source/'declared-inputs.npz');meta=json.loads((args.source/'protocol.json').read_text())
    weak=json.loads((args.source/'analysis.json').read_text())['weak_ids'];ids=sorted(set(j for i in weak for j in [2*(i//2),2*(i//2)+1]))
    jobs=[dict(id=i,parameters=d['parameters'][i].tolist(),initial=d['final_states'][i].tolist(),scale=inputs['scales'][i//2].tolist(),
        omega=float(inputs['omega'][i//2]),physical=meta['physical_parameters'][i//2],output=str(out)) for i in ids]
    write(out/'protocol.json',dict(source=str(args.source),sha256=hashlib.sha256((args.source/'complete-or-budgeted.npz').read_bytes()).hexdigest(),
        jobs=jobs,selection='Both starts at every tuple with a completed all-four-window rate>.001.',
        method='Independent literal physical field mapped to log coordinates,Radau2e-11/2e-13,maxstep.03,400xi.',
        return_rule='At least51 positive s=s_equilibrium crossings; maximum consecutive full3D log-coordinate return difference over last50<1e-7.',
        scope='Finite numerical period check,not a global exclusion or proof of physical validity. Mean-mass radius is not a bound on the lognormal size tail.'))
    for name in ['ice_cloud.py','check_ice_cloud_returns.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(run,j) for j in jobs]):
            row=f.result();rows.append(row);write(out/'analysis.json',dict(completed=len(rows),declared=len(jobs),runs=sorted(rows,key=lambda x:x['id'])))
            print(row,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

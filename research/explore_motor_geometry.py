"""Resolve geometry of retained motor chaos candidates without changing equations.

Fresh output directories only. Curves use accepted adaptive solver points;
their sampling density must not be interpreted as an invariant density.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from check_motor_pair_gpu_candidate import field as pair_field, physical
from motor_pair import original_field

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')


def fast_pair(x, p):
    q, y, z = x
    rho, r = p
    v = r*(q+y)
    return np.array([-rho*(q+y), v*(1-z)-y, v*y-z])


def pair_runs(args):
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    source = ROOT/'results/motor-pair/gpu-refined-complete.npz'
    z = np.load(source)
    ids = [437, 440, 441, 444]
    jobs = []
    for sid in ids:
        i = list(z['source_ids']).index(sid)
        jobs.append(dict(source_id=sid, parameters=z['parameters'][i].tolist(),
                         initial=z['balanced_states'][i].tolist(),
                         inherited_window_estimates=z['lle_windows'][i].tolist()))
    protocol = dict(source=source.relative_to(ROOT).as_posix(), source_sha256=sha(source),
                    script_sha256=sha(__file__), selection='One representative of each of the four retained r values; first source ID at each value.',
                    jobs=jobs, duration=args.duration, method='DOP853', rtol=1e-11,
                    atol=1e-13, max_step=args.step, chunks=100.,
                    sampling='All accepted solver points, no interpolation, no trajectory clipping.',
                    status='Numerical geometry follow-up; no new attractor, topology or priority claim.')
    write(out/'protocol.json', protocol)
    (out/'source.py').write_bytes(Path(__file__).read_bytes())
    records=[]
    for job in jobs:
        began=time.perf_counter(); p=np.array(job['parameters']); x=np.array(job['initial'])
        f,j=pair_field(x,p)
        error=float(np.max(abs(f-fast_pair(x,p))))
        transform=np.array([[1/p[1],-1,0],[0,1,0],[0,0,1.]])
        source_error=float(np.max(abs(f-transform@original_field(physical(x,p),p))))
        control_time=np.linspace(0,2,1001)
        controls={}
        for method in ['DOP853','Radau']:
            sol=solve_ivp(lambda t,s:fast_pair(s,p),(0,2),x,method=method,
                          t_eval=control_time,rtol=2e-12,atol=2e-14,max_step=.002)
            if not sol.success:raise RuntimeError(sol.message)
            controls[method]=sol.y.T
        control_error=float(np.max(abs(controls['DOP853']-controls['Radau'])))
        if max(error,source_error)>1e-10 or control_error>1e-7:raise RuntimeError('Control failed')
        times=[np.array([0.])]; states=[x[None]]; sections=[]; section_times=[]
        def section(t,s):return s[0]+s[1]
        section.direction=1
        elapsed=0.; count=0
        while elapsed<args.duration:
            end=min(elapsed+100,args.duration)
            sol=solve_ivp(lambda t,s:fast_pair(s,p),(elapsed,end),x,method='DOP853',
                          rtol=1e-11,atol=1e-13,max_step=args.step,events=section)
            if not sol.success:raise RuntimeError(sol.message)
            times.append(sol.t[1:]);states.append(sol.y[:,1:].T)
            sections.extend(sol.y_events[0]);section_times.extend(sol.t_events[0])
            x=sol.y[:,-1];elapsed=end;count+=sol.nfev
            print('pair',job['source_id'],'tau',elapsed,'seconds',round(time.perf_counter()-began,1),flush=True)
        t=np.concatenate(times); xs=np.vstack(states)
        path=out/f"pair-{job['source_id']}.npz"
        np.savez_compressed(path,time=t,states=xs,parameters=p,section_time=section_times,
                            section_states=np.array(sections).reshape(-1,3),control_time=control_time,
                            control_dop853=controls['DOP853'],control_radau=controls['Radau'])
        radius=xs[:,1]**2+(1-xs[:,2])**2
        assert np.isfinite(xs).all() and max(abs(xs[:,0]))<1+1e-8 and max(radius)<1+1e-8
        records.append(dict(source_id=job['source_id'],file=path.name,sha256=sha(path),
                            samples=len(t),sections=len(sections),nfev=count,
                            field_error=error,original_equation_error=source_error,
                            DOP853_Radau_max_state_error=control_error,
                            max_disk_radius_squared=float(max(radius)),max_abs_Q=float(max(abs(xs[:,0]))),
                            wall_seconds=time.perf_counter()-began))
        write(out/'analysis.json',dict(completed=len(records),declared=len(jobs),runs=records))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--duration',type=float,default=800.)
    parser.add_argument('--step',type=float,default=.005)
    args=parser.parse_args()
    pair_runs(args)

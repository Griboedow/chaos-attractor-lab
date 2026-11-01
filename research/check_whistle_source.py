"""Independent algebra and source-regime checks before any whistle search."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from whistle import reference,normalized,field_jac,dimensional_field,SOURCE


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);began=time.perf_counter()
    p=reference('B');coeff,scale,omega,derived=normalized(p)
    protocol=dict(source=SOURCE,equations='2,7-15;Table1-3',parameters=p,derived=derived,time='tau=omega_H*t',
        initial=[0.]*6,initial_pressure='pc=0 assumed; source specifies zero for the other five states but not pc explicitly.',
        ambient='rho=1.2041kg/m3,c=343.2m/s,nu=1.535e-5m2/s declared20C/1atm choices; exact numerical ambient constants not recovered.',
        source_access='Author-uploaded full text read on ResearchGate; publisher/institution PDF access unavailable. Printed equation signs retained.',
        observation_tau=args.duration,scope='Source implementation/control only, no broad scan or chaos claim.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['whistle.py','check_whistle_source.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rng=np.random.default_rng(2026110105);algebra=0.;jacerror=0.
    for _ in range(30):
        q=rng.uniform(-.7,.7,6);q[[0,1,2]]+=1.
        f,j=field_jac(q,coeff);independent=dimensional_field(q*scale,p)/(scale*omega)
        algebra=max(algebra,float(np.max(abs(f-independent)/(1+abs(independent)))))
        cols=[]
        for k in range(6):
            delta=1e-5*np.eye(6)[k]
            cols.append((dimensional_field((q+delta)*scale,p)-dimensional_field((q-delta)*scale,p))/(scale*omega*2e-5))
        finite=np.column_stack(cols);jacerror=max(jacerror,float(np.max(abs(j-finite)/(1+abs(j)))))
    assert algebra<1e-11 and jacerror<1e-6
    q=np.zeros(6);grid=np.linspace(0,100,5001)
    explicit=solve_ivp(lambda t,q:field_jac(q,coeff)[0],(0,100),q,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.02,t_eval=grid)
    implicit=solve_ivp(lambda t,q:dimensional_field(q*scale,p)/(scale*omega),(0,100),q,
        method='Radau',rtol=2e-11,atol=2e-13,max_step=.02,t_eval=grid)
    assert explicit.success and implicit.success
    difference=float(np.max(abs(explicit.y-implicit.y)/(1+abs(implicit.y))));assert difference<1e-7
    np.savez_compressed(out/'independent-short-control.npz',time=grid,explicit=explicit.y.T,implicit=implicit.y.T,scale=scale,omega=omega,coefficients=coeff)
    print('controls',algebra,jacerror,difference,derived,flush=True)
    # Long source-zero-start history; the early comparison is an independent
    # implementation control, not a substitute for long solver convergence.
    times=np.linspace(0,args.duration,1+int(round(args.duration/.05)))
    sol=solve_ivp(lambda t,q:field_jac(q,coeff)[0],(0,args.duration),q,method='Radau',jac=lambda t,q:field_jac(q,coeff)[1],
        rtol=2e-9,atol=2e-12,max_step=.1,t_eval=times)
    assert sol.success
    path=out/'source-history.npz';np.savez_compressed(path,time=times,states=sol.y.T,scale=scale,omega=omega,coefficients=coeff)
    tail=sol.y[:,times>args.duration*.9];pressure=tail[2]*scale[2];vortex=tail[4]
    record=dict(algebra_error=algebra,jacobian_error=jacerror,independent_short_solver_error=difference,
        duration_tau=args.duration,duration_seconds=args.duration/omega,tail_pressure_mean=float(np.mean(pressure)),
        tail_pressure_peak_to_peak=float(np.ptp(pressure)),tail_vortex_peak_to_peak=float(np.ptp(vortex)),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),wall_seconds=time.perf_counter()-began,
        limitation='Long floating-point source-zero-start history with declared ambient/initial-pressure assumptions; not yet a calibrated experimental reproduction.')
    (out/'analysis.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--duration',type=float,default=30000.);main(p.parse_args())

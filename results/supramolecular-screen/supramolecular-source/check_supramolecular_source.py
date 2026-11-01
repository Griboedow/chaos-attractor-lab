"""Source-cycle and literal four-state checks before a polymer search."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from supramolecular import REFERENCE,field_jac,independent_field,source_four_field,equilibrium


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);p=REFERENCE.copy();eq=equilibrium(p)
    initial=eq['state']*np.array([1.01,.999,1.005]);beta=p[1]
    original=np.array([beta*initial[1],beta*initial[0],beta*(p[2]-initial[0]-initial[1]),initial[2]])
    grid=np.linspace(0,100,20001)
    def section(t,x):return x[0]-1
    section.direction=1
    dop=solve_ivp(lambda t,x:field_jac(x,p)[0],(0,100),initial,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.005,t_eval=grid,events=section)
    rad=solve_ivp(lambda t,v:source_four_field(v,p),(0,100),original,method='Radau',rtol=2e-11,atol=2e-13,max_step=.005,t_eval=grid)
    assert dop.success and rad.success
    transformed=np.array([rad.y[1]/beta,rad.y[0]/beta,rad.y[3]])
    agreement=float(np.max(abs(dop.y-transformed)/(1+abs(transformed))))
    balance=float(max(abs(rad.y[0]+rad.y[1]+rad.y[2]-beta*p[2])))
    algebra=0.
    for x in dop.y.T[::1000]:
        f,j=field_jac(x,p);jj=np.column_stack([independent_field(x+1e-25j*np.eye(3)[k],p).imag/1e-25 for k in range(3)])
        algebra=max(algebra,float(np.max(abs(j-jj)/(1+abs(jj)))),float(np.max(abs(f-independent_field(x,p))/(1+abs(f)))))
    assert agreement<1e-6 and algebra<1e-10 and balance<1e-6
    events=dop.t_events[0];tail=events[events>80];period=float(np.median(np.diff(tail)));assert .7<period<.8
    x=dop.y.T;length=beta*(p[2]-x[:,0]-x[:,1])/x[:,2]
    np.savez_compressed(out/'source-control.npz',parameters=p,time=grid,states=x,source_states=rad.y.T,initial=initial,section_times=events,section_states=dop.y_events[0])
    report=dict(source='Sharko2022 simplified moment closure; thesis pp48-49 Eqs5.6/5.12 and Figure14; k+=kd3=1,ctot10000,ka.2,kd21000,C4.6.',
        parameters=p.tolist(),equilibrium=eq['state'].tolist(),equilibrium_eigenvalues=[[float(v.real),float(v.imag)] for v in eq['eigenvalues']],
        independent_four_state_solver_scaled_error=agreement,source_mass_balance_error=balance,source_jacobian_relative_error=algebra,
        period_tau=period,source_rounded_period=.75,observed_length_range=[float(length.min()),float(length.max())],
        warning='Approximate moment model for long polymers; not the full size-distribution equations. No chaos or novelty claim.',
        sha256=hashlib.sha256((out/'source-control.npz').read_bytes()).hexdigest())
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    for name in ['supramolecular.py','check_supramolecular_source.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

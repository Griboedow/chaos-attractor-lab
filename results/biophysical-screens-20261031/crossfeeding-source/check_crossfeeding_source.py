"""Source chemostat reproduction, independent solver and budget controls."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from crossfeeding import REFERENCE,REFERENCE_START,SOURCE,AUTHOR_COMMIT,PARAMETERS,field_jac,independent_field,balance_weights,balance_in,branch_gaps


def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    p=REFERENCE;x=REFERENCE_START
    write(out/'protocol.json',dict(source=SOURCE,author_commit=AUTHOR_COMMIT,parameters=PARAMETERS,
        tuple=p.tolist(),initial=x.tolist(),time='tau=D*t; t in author hours; source time0..200 hours is tau0..40.',
        method='DOP853 versus independent stoichiometric Radau, rtol2e-11 atol2e-13 maxstep.005tau.',
        limitation='Continuous min branches; no smoothing, forcing, clipping or fitted parameter changes.'))
    for name in ['crossfeeding.py','check_crossfeeding_source.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    grid=np.linspace(0,40,20001)
    dop=solve_ivp(lambda t,y:field_jac(y,p)[0],(0,40),x,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.005,t_eval=grid)
    rad=solve_ivp(lambda t,y:independent_field(y,p),(0,40),x,method='Radau',rtol=2e-11,atol=2e-13,max_step=.005,t_eval=grid)
    assert dop.success and rad.success
    states=dop.y.T;error=float(np.max(abs(dop.y-rad.y)/(1+abs(rad.y))))
    algebra=0.;budget_field=0.;weights=balance_weights(p);supply=balance_in(p)
    for state in states[::997]:
        assert min(abs(branch_gaps(state,p)))>1e-9
        f,j=field_jac(state,p)
        jj=np.column_stack([independent_field(state+1e-25j*np.eye(5)[i],p).imag/1e-25 for i in range(5)])
        algebra=max(algebra,float(np.max(abs(j-jj)/(1+abs(jj)))),float(np.max(abs(f-independent_field(state,p)))))
        budget_field=max(budget_field,float(abs(weights@f-(supply-weights@state))))
    expected=supply+(weights@x-supply)*np.exp(-grid)
    budget_error=float(np.max(abs(states@weights-expected)))
    assert error<2e-6 and algebra<1e-9 and budget_field<1e-9 and budget_error<1e-8
    path=out/'source-control.npz';np.savez_compressed(path,time_tau=grid,time_hours=grid/p[12],states=states,parameters=p,initial=x)
    record=dict(independent_solver_scaled_max_error=error,algebra_error=algebra,
        budget_field_residual=budget_field,budget_history_max_residual=budget_error,budget_weights=weights.tolist(),budget_in=supply,
        final10tau_minimum=states[grid>=30].min(axis=0).tolist(),final10tau_maximum=states[grid>=30].max(axis=0).tolist(),
        samples=len(grid),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    write(out/'analysis.json',record);print(json.dumps(record),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

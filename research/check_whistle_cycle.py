"""Numerical period-one control at the literal whistle source tuple."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from whistle import reference,normalized,field_jac,dimensional_field


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);p=reference('B');coeff,scale,omega,derived=normalized(p)
    assert np.allclose(coeff,source['coefficients'],rtol=0,atol=0)
    def event(t,q):return q[4]
    event.direction=1
    sol=solve_ivp(lambda t,q:field_jac(q,coeff)[0],(0,100),source['states'][-1],method='Radau',
        jac=lambda t,q:field_jac(q,coeff)[1],rtol=2e-11,atol=2e-13,max_step=.02,events=event)
    assert sol.success and len(sol.t_events[0])>=5
    initial=sol.y_events[0][-1];period=np.diff(sol.t_events[0])[-1];indices=[0,1,2,3,5]
    def shoot(z,return_all=False):
        q=np.zeros(6);q[indices]=z[:-1];T=z[-1]
        def augmented(t,y):
            f,j=field_jac(y[:6],coeff)
            return np.r_[f,(j@y[6:].reshape(6,6)).ravel()]
        s=solve_ivp(augmented,(0,T),np.r_[q,np.eye(6).ravel()],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.005)
        assert s.success
        end=s.y[:6,-1];M=s.y[6:,-1].reshape(6,6);f=field_jac(end,coeff)[0]
        if return_all:return q,T,end,M
        return end-q,np.column_stack([(M-np.eye(6))[:,indices],f])
    result=root(lambda z:shoot(z),np.r_[initial[indices],period],jac=True,tol=1e-10)
    q,T,end,M=shoot(result.x,True);error=float(np.max(abs(end-q)))
    assert result.success and error<1e-8 and T>0
    grid=np.linspace(0,T,4001)
    independent=solve_ivp(lambda t,x:dimensional_field(x*scale,p)/(scale*omega),(0,T),q,method='Radau',
        rtol=2e-12,atol=2e-14,max_step=.002,t_eval=grid)
    assert independent.success
    vals=np.linalg.eigvals(M);phase=int(np.argmin(abs(vals-1)));other=np.delete(vals,phase)
    history=independent.y.T;pressure=history[:-1,2]*scale[2];ac=pressure-pressure.mean()
    record=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        parameters=p,derived=derived,period_tau=float(T),frequency_Hz=float(omega/T),observed_source_Hz=p['observed_frequency'],
        relative_frequency_difference=float(omega/T/p['observed_frequency']-1),shooting_residual=error,
        independent_dimensional_Radau_closure=float(np.max(abs(independent.y[:,-1]-q))),
        phase_multiplier=[float(vals[phase].real),float(vals[phase].imag)],
        nonflow_multipliers=[[float(v.real),float(v.imag)] for v in other],
        largest_nonflow_magnitude=float(np.max(abs(other))),pressure_mean_Pa=float(pressure.mean()),
        pressure_peak_to_peak_Pa=float(np.ptp(pressure)),pressure_AC_rms_Pa=float(np.sqrt(np.mean(ac**2))),
        state_extrema=np.array([history.min(0),history.max(0)]).tolist(),
        interpretation='Numerically attracting source cycle in literal equations with declared ambient/pc-initial assumptions; no chaos or exact experimental calibration claim. Tiny strongly contracting multipliers are not relatively accurate.',
        frequency_limit='Empirical Strouhal frequency is fitted in the source; matching the experimental peak is not independent validation of all dynamics.')
    (out/'analysis.json').write_text(json.dumps(record,indent=2)+'\n')
    np.savez_compressed(out/'cycle.npz',time=grid,states=history,initial=q,period=T,monodromy=M,coefficients=coeff,scale=scale,omega=omega)
    for name in ['check_whistle_cycle.py','whistle.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

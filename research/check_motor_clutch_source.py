"""Independent source-form controls before any motor-clutch parameter survey."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from motor_clutch import SOURCE,NAMES,parameters,field_jac,independent_field,equilibria


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    cases=[parameters(vu=v) for v in [15.,3.,1.5,.5]]
    report=[]
    (out/'protocol.json').write_text(json.dumps(dict(source=SOURCE,equations='Literal smooth Eq6/A1',
        parameters=NAMES,tuples=[p.tolist() for p in cases],
        comparison='Figure3a: ks=kc, gamma=Gamma. A1 does not list precise vu values; these are declared representative controls, not exact curve reproductions.',
        force_scale='Numerical TableI f=2.3033 interpreted as pN; its printed unit is inconsistent with the defining scale.',
        limitations='Eq6 continues the linear motor-velocity branch beyond stall; source physical piecewise law is a separate interpretation.'),indent=2)+'\n')
    for name in ['motor_clutch.py','check_motor_clutch_source.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for i,p in enumerate(cases):
        eq=equilibria(p)[-1];x=eq['state']*np.array([1.01,1.002,1.001,1.005,1.001])
        f,j=field_jac(x,p)
        jj=np.column_stack([independent_field(x+1e-25j*np.eye(5)[k],p).imag/1e-25 for k in range(5)])
        algebra=max(np.max(abs(f-independent_field(x,p))),np.max(abs(j-jj)))
        assert algebra<1e-10
        time=np.linspace(0,4000,40001)
        dop=solve_ivp(lambda t,y:field_jac(y,p)[0],(0,4000),x,method='DOP853',rtol=2e-10,atol=2e-12,max_step=.2,t_eval=time)
        rad=solve_ivp(lambda t,y:independent_field(y,p),(0,4000),x,method='Radau',rtol=2e-10,atol=2e-12,max_step=1.,t_eval=time)
        assert dop.success and rad.success
        error=float(np.max(abs(dop.y-rad.y)/(1+abs(rad.y))))
        assert error<1e-5
        states=dop.y.T;tail=states[30000:]
        path=out/f'case-{i}.npz';np.savez_compressed(path,time=time,states=states,parameters=p,initial=x)
        record=dict(id=i,parameters=p.tolist(),equilibrium=eq['state'].tolist(),
            eigenvalues=[[float(z.real),float(z.imag)] for z in eq['eigenvalues']],
            equilibrium_residual=eq['residual'],algebra_error=float(algebra),independent_state_error=error,
            tail_minimum=tail.min(axis=0).tolist(),tail_maximum=tail.max(axis=0).tolist(),
            superstall_sample_fraction=float(np.mean(states[:,3]*p[5]>p[8])),
            negative_clutch_extension_fraction=float(np.mean(states[:,2]<states[:,4])),
            file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        report.append(record);(out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
        print(i,'vu',p[11],'agreement',error,'span',np.ptp(tail,axis=0),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())

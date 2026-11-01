"""Independent source algebra and solver controls; no fitted-amplitude claim."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from wing_dsa import field_jac,source_equivalence,equilibrium_polynomial,SOURCE


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(2026110102);algebra=0.;jac_error=0.;poly_error=0.
    for i in range(40):
        p=np.array([10**rng.uniform(-1,1),rng.uniform(.01,.99),10**rng.uniform(-2,2),1.])
        q=rng.normal(size=4);f,j=field_jac(q,p)
        algebra=max(algebra,float(np.max(abs(f-source_equivalence(q,p))/(1+abs(f)))))
        # C1 quadratic drag has derivative 2|v|; use centered differences away from v=0.
        eps=1e-5
        jj=np.column_stack([(source_equivalence(q+eps*np.eye(4)[n],p)-source_equivalence(q-eps*np.eye(4)[n],p))/(2*eps) for n in range(4)])
        jac_error=max(jac_error,float(np.max(abs(j-jj)/(1+abs(j)))))
        actual=np.poly(field_jac(np.zeros(4),p)[1]);expected=equilibrium_polynomial(p)
        poly_error=max(poly_error,float(np.max(abs(actual-expected)/(1+abs(expected)))))
    assert algebra<1e-12 and jac_error<1e-8 and poly_error<1e-11
    controls=[('stable',[.5,.5,.2,1.]),('beyond_threshold',[.5,.5,4.,1.]),
        ('reported_mechanics',[35/np.sqrt(.416/.00196),.5,.06/.416,1.])]
    records=[]
    for name,parameters in controls:
        p=np.array(parameters);initial=np.array([.1,0,0,0.]);rhs=lambda t,x:field_jac(x,p)[0]
        grid=np.linspace(0,1600,80001)
        one=solve_ivp(rhs,(0,1600),initial,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.04,t_eval=grid)
        two=solve_ivp(lambda t,x:source_equivalence(x,p),(0,1600),initial,method='Radau',
            jac=lambda t,x:field_jac(x,p)[1],rtol=2e-10,atol=2e-12,max_step=.04,t_eval=grid)
        assert one.success and two.success
        difference=float(np.max(abs(one.y-two.y)/(1+abs(one.y))))
        assert difference<2e-5
        path=out/f'{name}.npz';np.savez_compressed(path,time=grid,states=one.y.T,independent_states=two.y.T,parameters=p,initial=initial)
        records.append(dict(name=name,parameters=p.tolist(),solver_scaled_difference=difference,
            equilibrium_eigenvalues=[[float(z.real),float(z.imag)] for z in np.linalg.eigvals(field_jac(np.zeros(4),p)[1])],
            tail_amplitude=float(np.ptp(one.y[0,grid>=1200])/2),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        print(name,difference,records[-1]['tail_amplitude'],flush=True)
    report=dict(source=SOURCE,equations='4,9 and12',controls=records,algebra_error=algebra,jacobian_error=jac_error,characteristic_polynomial_error=poly_error,
        time='tau=sqrt(k/I)*t',scope='Exact dimensionless mathematical model. Quadratic drag sets the amplitude unit and is not inferred from a source figure.',
        physical_control='k=.416Nm/rad,I=.00196kg m^2,r3=35/s,kappa=.5,mu=.06 as reported mechanical/gain values. Drag=.003 is an arbitrary algebra-test unit, not a fitted source value.',
        limitations='No quantitative experimental-amplitude reproduction. The 2022 paper neglects bearing friction; printed Eq10 omits drag division by inertia, so literal dimensional Eq9 is used.')
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    for name in ['wing_dsa.py','check_wing_dsa_source.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

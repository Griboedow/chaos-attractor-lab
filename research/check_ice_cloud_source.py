"""Independent dimensional controls before any ice-cloud chaos screen."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from ice_cloud import coefficients,make_normalized,field_jac,physical_field,to_physical,SOURCE


def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ['ice_cloud.py','check_ice_cloud_source.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    write(out/'protocol.json',dict(source=SOURCE,model='Full Eq55--57 including nonzero nucleation mass and piecewise evaporation; no n_small regularization.',
        controls='Three published low/middle/high updraft examples at declared T230K,p300hPa,F1,Na300/cm3,cp1005J/kg/K. Reproduces qualitative branch types, not exact source figure traces.',
        normalization='u=log(n/n0),v=log(q/q0),s; full positive equilibrium n0,q0; tau=As*t.',
        caveats=['Specific heat cp1005 is a declared standard convention; no numerical value or author implementation was recovered.',
                 'The source sometimes uses Na10000/cm3 for reference comparisons. Here Na300/cm3 is explicit.',
                 'The final2026 article was identified, but this implementation is pinned to the accessible2025 preprint v1.',
                 'Constant p,T with updraft forcing is the published autonomous approximation, not a rising parcel thermodynamic history.']))
    rng=np.random.default_rng(2026110111);algebra=[];jacerr=[];controls=[]
    for w in [.01,.1,1.]:
        c=coefficients(updraft=w);p,scale,omega=make_normalized(c)
        for _ in range(16):
            z=np.r_[rng.uniform(-2,2,2),rng.uniform(.85,1.6)]
            x=to_physical(z,scale);f,j=field_jac(z,p)
            direct=physical_field(x,c)/np.r_[x[:2],1.]/omega
            algebra.append(float(np.max(abs(f-direct)/(1+abs(direct)))))
            eps=1e-5;eye=np.eye(3)*eps
            finite=np.column_stack([(-field_jac(z+2*h,p)[0]+8*field_jac(z+h,p)[0]-8*field_jac(z-h,p)[0]+field_jac(z-2*h,p)[0])/(12*eps) for h in eye])
            jacerr.append(float(np.max(abs(j-finite)/(1+abs(j)))))
        z0=np.array([np.log(.99),np.log(.99),p[-1]*.99]);x0=to_physical(z0,scale)
        # Raw physical variables, with constant scales only, independent of the logarithmic implementation.
        raw0=x0/scale
        def raw(t,y):return physical_field(y*scale,c)/(scale*omega)
        short=solve_ivp(raw,(0,1),raw0,method='Radau',rtol=2e-11,atol=2e-13,max_step=.001,dense_output=True)
        ref=solve_ivp(lambda t,z:field_jac(z,p)[0],(0,1),z0,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.001,dense_output=True)
        assert short.success and ref.success
        times=np.linspace(0,1,1001);a=short.sol(times).T*scale;b=to_physical(ref.sol(times).T,scale)
        error=float(np.max(abs((a-b)/scale)/(1+abs(a/scale))))
        assert error<1e-6
        def section(t,z):return z[2]-p[-1]
        section.direction=1
        total=solve_ivp(lambda t,z:field_jac(z,p)[0],(0,100),z0,method='Radau',jac=lambda t,z:field_jac(z,p)[1],
            rtol=2e-9,atol=2e-11,max_step=.01,events=section,dense_output=True)
        assert total.success
        sample_t=np.linspace(90,100,10001);z=total.sol(sample_t).T;sections=total.y_events[0];events=total.t_events[0]
        tail=sections[-101:];distinct=bool(np.ptp(z[:,2])>1e-6)
        recurrence=None
        if len(tail)>2 and distinct:
            recurrence=dict(last100_section_lag1_max=float(np.max(np.linalg.norm(np.diff(tail,axis=0),axis=1))),
                period_tau=float(np.mean(np.diff(events[-101:]))),period_seconds=float(np.mean(np.diff(events[-101:]))/omega))
        path=out/f'updraft-{w:g}.npz'
        np.savez_compressed(path,physical_coefficients=c,normalized_parameters=p,scale=scale,As=omega,
            time=sample_t,log_states=z,physical_states=to_physical(z,scale),section_times=events,section_states=sections,
            short_time=times,short_dimensional=a,short_log_coordinates=b)
        record=dict(updraft=w,As=omega,equilibrium=[scale[0],scale[1],p[-1]],
            equilibrium_eigenvalues_per_tau=[[float(v.real),float(v.imag)] for v in np.linalg.eigvals(field_jac([0,0,p[-1]],p)[1])],
            short_solver_scaled_error=error,tail_saturation_peak_to_peak=float(np.ptp(z[:,2])),
            tail_log_state_peak_to_peak=np.ptp(z,axis=0).tolist(),recurrence=recurrence,
            final_field_norm=float(np.linalg.norm(field_jac(z[-1],p)[0])),file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        controls.append(record);write(out/'progress.json',controls);print(record,flush=True)
    assert max(algebra)<1e-10 and max(jacerr)<1e-5
    write(out/'validation.json',dict(passed=True,max_scaled_algebra_error=max(algebra),max_scaled_five_point_jacobian_error=max(jacerr),runs=controls))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

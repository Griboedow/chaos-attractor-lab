"""Declared empirical-coefficient screen at fixed source-B geometry and flow.

No fitted Strouhal-law extrapolation: geometry, supply pressure and jet speed
are held at the published B example. Coefficient variations are exploratory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import qmc
from whistle import reference,normalized,field_jac,dimensional_field,SOURCE
from whistle_gpu import WhistleGPU


def write(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def save(out,result):
    arrays={k:v for k,v in result.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(out/'complete-or-budgeted.npz',**arrays)
    meta={k:v for k,v in result.items() if k not in arrays}
    write(out/'settings.json',meta)


def controls(out,gpu):
    source=reference('B');variants=[source.copy() for _ in range(5)]
    for p,cv,tau,gamma in zip(variants[1:],[.25,4.,4.,.25],[.25,4.,.25,4.],[.5,2.,.5,2.]):
        p['Cv']*=cv;p['tau0']*=tau;p['gamma']*=gamma
    rng=np.random.default_rng(2026110108);states=[];parameters=[];tangents=[]
    for p in variants:
        coeff,_,_,_=normalized(p)
        for _ in range(8):
            q=rng.uniform(-.5,.5,6);q[:3]+=.7
            states.append(q);parameters.append(coeff);tangents.append(rng.normal(size=6))
    states=np.array(states);parameters=np.array(parameters);tangents=np.array(tangents)
    actual=gpu.values(states,parameters,tangents);expected=[]
    for q,p,v in zip(states,parameters,tangents):
        f,j=field_jac(q,p);expected.append(np.r_[f,j@v])
    expected=np.array(expected)
    field_error=float(np.max(abs(actual-expected)/(1+abs(expected))));assert field_error<1e-11
    rows=[];arrays={}
    for i,p in enumerate(variants):
        coeff,scale,omega,_=normalized(p)
        x=np.sqrt(coeff[0]/(2*coeff[2]));q=np.array([x,x,.5,0.,2.,0.])
        tangent=np.sin(np.arange(6)+1.);tangent/=np.linalg.norm(tangent)
        def augmented(t,z):
            f=dimensional_field(z[:6]*scale,p)/(scale*omega);j=field_jac(z[:6],coeff)[1]
            return np.r_[f,j@z[6:]]
        sol=solve_ivp(augmented,(0,20),np.r_[q,tangent],method='Radau',rtol=2e-11,atol=2e-13,max_step=.01)
        assert sol.success
        v=sol.y[6:,-1].copy();rate=float(np.log(np.linalg.norm(v))/20);v/=np.linalg.norm(v)
        # Also compare a second independent integration method.
        sol2=solve_ivp(augmented,(0,20),np.r_[q,tangent],method='DOP853',rtol=2e-11,atol=2e-13,max_step=.005)
        assert sol2.success
        ref=sol2.y[6:,-1].copy();rate2=float(np.log(np.linalg.norm(ref))/20);ref/=np.linalg.norm(ref)
        cpu_difference=float(np.max(abs(sol.y[:6,-1]-sol2.y[:6,-1])/(1+abs(sol.y[:6,-1]))))
        for tol in [1e-7,1e-9]:
            result=gpu.scan(coeff[None],q[None],warm=0,duration=20,rtol=tol,max_step=.03,max_wall=120)
            row=dict(control=i,physical_parameters=p,rtol=tol,status=int(result['status'][0]),
                cpu_state_difference=cpu_difference,state_error=float(np.max(abs(result['final_states'][0]-sol.y[:6,-1])/(1+abs(sol.y[:6,-1])))),
                tangent_error=float(np.max(abs(result['final_tangents'][0]-v))),
                cpu_tangent_direction_difference=float(np.max(abs(ref-v))),cpu_rate_difference=abs(rate-rate2),cpu_rate=rate,gpu_rate=float(result['lle'][0]),
                rate_error=float(abs(rate-result['lle'][0])),counts=result['step_counts'][0].tolist())
            rows.append(row);print({k:v for k,v in row.items() if k!='physical_parameters'},flush=True)
            assert row['status']==0 and row['cpu_state_difference']<1e-7 and row['state_error']<1e-5 and row['tangent_error']<1e-4 and row['rate_error']<1e-5
        arrays[f'control_{i}']=np.r_[q,sol.y[:6,-1],v]
    write(out/'validation.json',dict(field_tangent_scaled_error=field_error,rows=rows,passed=True,
        scope='Independent dimensional state implementation, Radau and DOP853; piecewise analytic tangent. Short controls do not classify long dynamics.'))
    np.savez_compressed(out/'controls.npz',**arrays)


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ['search_whistle.py','whistle.py','whistle_gpu.py','whistle.cl']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    gpu=WhistleGPU()
    if args.validate_only:controls(out,gpu);return
    validated=json.loads(args.validation.read_text());assert validated['passed']
    columns=['Cv','tau0','gamma','beta','mu','kappa'];lower=np.array([.25,.25,.5,.5,1.,.25]);upper=np.array([4.,4.,2.,2.,4.,4.])
    points=qmc.Sobol(6,scramble=True,seed=2026110107).random_base2(8)
    multipliers=np.vstack([np.ones(6),np.exp(np.log(lower)+points*np.log(upper/lower))])
    physical=[];coefficients=[];initial=[]
    for row in multipliers:
        p=reference('B')
        for name,factor in zip(columns,row):p[name]*=float(factor)
        c,_,_,_=normalized(p);x=(-c[1]+np.sqrt(c[1]**2+2*c[2]*c[0]))/(2*c[2])
        physical.append(p)
        for phase in range(2):
            coefficients.append(c)
            initial.append([x,x,.5,0.,2. if phase==0 else 0.,0. if phase==0 else 2.*np.sqrt(c[9])])
    coefficients=np.array(coefficients);initial=np.array(initial)
    protocol=dict(source=SOURCE,validation=str(args.validation),validation_sha256=hashlib.sha256(args.validation.read_bytes()).hexdigest(),
        model='Literal six-state Eq2,7,13,15 with empirical Table3 coefficients; fixed B geometry,pressure93Pa,velocity8.16m/s.',
        columns=columns,multiplier_lower=lower.tolist(),multiplier_upper=upper.tolist(),physical_parameters=physical,
        tuples=len(multipliers),lanes=len(initial),seed=2026110107,
        starts='Two prescribed vortex phases at amplitude2,acoustic dc equilibrium,g=0. Not starts on a numerically fitted attractor.',
        warm=1000.,duration=4000.,rtol=1e-8,max_step=.03,
        primary_rule='Completed and all four1000tau leading estimates>.002.',weak_rule='Completed and all four leading estimates>.0005.',
        physical_scope='Empirical coefficient sensitivity, not experimentally calibrated new apparatus or validity across all coefficient values. Fitted St is not extrapolated to new flow/geometry.',
        novelty='No exact-model or aeroacoustic-chaos priority claim; existing van der Pol wake-oscillator equivalents require comparison.')
    write(out/'protocol.json',protocol)
    np.savez_compressed(out/'declared-inputs.npz',parameters=coefficients,initial=initial,multipliers=multipliers,
        tuple_id=np.repeat(np.arange(len(multipliers)),2),start_id=np.tile([0,1],len(multipliers)))
    result=gpu.scan(coefficients,initial)
    result.update(tuple_id=np.repeat(np.arange(len(multipliers)),2),start_id=np.tile([0,1],len(multipliers)))
    save(out,result);status=result['status'];good=status==0
    retained=np.flatnonzero(good&(result['lle_windows'].min(1)>.002));weak=np.flatnonzero(good&(result['lle_windows'].min(1)>.0005))
    write(out/'analysis.json',dict(lanes=len(initial),statuses={str(int(k)):int(v) for k,v in zip(*np.unique(status,return_counts=True))},
        retained_ids=retained.tolist(),weak_ids=weak.tolist(),wall_seconds=result['wall_seconds'],
        limitation='Finite exploratory screen; failures unresolved; positive rates require independent full-spectrum and attraction checks.'))
    print('retained',retained,'weak',weak,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--validate-only',action='store_true');p.add_argument('--validation',type=Path);main(p.parse_args())

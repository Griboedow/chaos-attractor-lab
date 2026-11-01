"""Ice-cloud coefficient-family screen, retaining the full nucleation terms."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import qmc
from ice_cloud import coefficients,make_normalized,field_jac,physical_field,to_physical,SOURCE
from ice_cloud_gpu import IceCloudGPU


def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')


def prepare(physical):
    c=coefficients(**physical);p,scale,As=make_normalized(c)
    eig=np.linalg.eigvals(field_jac(np.array([0.,0.,p[-1]]),p)[1]);imag=max(abs(eig.imag))
    period=2*np.pi/(imag if imag>1e-6 else max(abs(eig)))
    p=p.copy();p[:8]*=period
    return c,p,scale,As/period,period,eig


def controls(out,gpu):
    variants=[dict(T=230.,pressure=30000.,updraft=w,sedimentation=1.,aerosols=3e8,r0=3.) for w in [.01,.1,1.]]
    variants += [dict(T=200.,pressure=25000.,updraft=.01,sedimentation=.1,aerosols=3e7,r0=1.2),
                 dict(T=230.,pressure=22000.,updraft=.5,sedimentation=.1,aerosols=3e9,r0=20.)]
    states=[];params=[];vectors=[];rng=np.random.default_rng(2026110113)
    for physical in variants:
        c,p,scale,omega,_,_=prepare(physical)
        for _ in range(10):
            states.append(np.r_[rng.uniform(-2.,2.,2),rng.uniform(.9,1.6)])
            params.append(p);vectors.append(rng.normal(size=3))
    actual=gpu.values(np.array(states),np.array(params),np.array(vectors));expected=[]
    for z,p,v in zip(states,params,vectors):
        f,j=field_jac(z,p);expected.append(np.r_[f,j@v])
    expected=np.array(expected);err=float(np.max(abs(actual-expected)/(1+abs(expected))));assert err<1e-10
    records=[]
    for physical in variants:
        c,p,scale,omega,_,_=prepare(physical)
        z0=np.array([np.log(.99),np.log(.99),p[-1]*.99]);v=np.sin(np.arange(3)+1.);v/=np.linalg.norm(v)
        def augmented(t,z):
            q=to_physical(z[:3],scale)
            f=physical_field(q,c)/(np.r_[q[:2],1.]*omega)
            return np.r_[f,field_jac(z[:3],p)[1]@z[3:]]
        sol=solve_ivp(augmented,(0,2),np.r_[z0,v],method='Radau',rtol=2e-11,atol=2e-13,max_step=.002)
        assert sol.success
        norm=np.linalg.norm(sol.y[3:,-1]);direction=sol.y[3:,-1]/norm;rate=np.log(norm)/2
        for tol in [1e-8,1e-10]:
            r=gpu.scan(p[None],z0[None],warm=0,duration=2,rtol=tol,max_step=.01,max_wall=60)
            record=dict(physical=physical,rtol=tol,status=int(r['status'][0]),
                state_error=float(np.max(abs(r['final_states'][0]-sol.y[:3,-1]))),
                tangent_error=float(np.max(abs(r['final_tangents'][0]-direction))),rate_error=float(abs(r['lle'][0]-rate)))
            assert record['status']==0 and max(record[k] for k in ['state_error','tangent_error','rate_error'])<2e-6
            records.append(record);print(record,flush=True)
    write(out/'validation.json',dict(passed=True,field_tangent_scaled_error=err,records=records,
        scope='Independent dimensional physical field,Radau,analytic tangent,FP64 GPU controls; finite short intervals only.'))


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ['ice_cloud.py','ice_cloud_gpu.py','ice_cloud.cl','search_ice_cloud.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    gpu=IceCloudGPU()
    if args.validate_only:controls(out,gpu);return
    validated=json.loads(args.validation.read_text());assert validated['passed']
    # Vary the previously fixed mass-distribution shape and aerosol abundance,
    # alongside the original environmental control family. No new state or law.
    points=qmc.Sobol(6,scramble=True,seed=2026110112).random_base2(10)
    physical=[];normalized=[];initial=[];scales=[];omegas=[];equilibria=[]
    for a in points:
        phys=dict(T=190+50*a[0],pressure=20000+10000*a[1],updraft=.0005*4000**a[2],
            sedimentation=.01*100**a[3],aerosols=3e6*1000**a[4],r0=1.05*(30/1.05)**a[5])
        c,p,scale,omega,period,eig=prepare(phys)
        physical.append(phys);scales.append(scale);omegas.append(omega);equilibria.append([scale[0],scale[1],p[-1]])
        normalized.extend([p,p]);initial.extend([[np.log(.99),np.log(.99),p[-1]*.99],[np.log(.2),np.log(5.),1.2]])
    normalized=np.array(normalized);initial=np.array(initial)
    protocol=dict(source=SOURCE,validation=str(args.validation),validation_sha256=hashlib.sha256(args.validation.read_bytes()).hexdigest(),
        physical_parameters=physical,tuples=len(physical),lanes=len(initial),seed=2026110112,
        source_family='Full preprint Eq55--57; r0 and aerosol abundance now varied through the printed coefficient formulas. No artificial forcing or added state.',
        ranges=dict(T_K=[190,240],pressure_Pa=[20000,30000],updraft_m_s=[.0005,2],sedimentation_F=[.01,1],aerosols_m3=[3e6,3e9],r0=[1.05,30]),
        convention='cp1005J/kg/K; remaining constants exactly as ice_cloud.py. Source distribution width r0=3; broad width variation is exploratory and needs physical validity checks.',
        starts='Near full equilibrium (all physical components*0.99); separate n*.2,q*5,s1.2.',
        time='xi=As*t/Tlin; Tlin=2*pi/imaginary eigenfrequency in tau=As*t, or2*pi/max(abs(eigenvalues)) on all-real branches. Tuple-dependent fixed positive scaling only.',
        warm=200.,duration=800.,rtol=1e-9,max_step=.03,primary='Completed with all four200xi leading estimates>.005.',weak='Completed with all four estimates>.001.',
        validity='Numerical log cutoff abs650 or saturation outside0..5 is unresolved. Must assess negligible aerosol depletion,crystal size,parameter validity and solver convergence before interpreting a positive signal.',
        novelty='No absence-of-prior-chaos or new-topology claim. The original study already maps periodic cycles in the fixed-width environmental family.')
    write(out/'protocol.json',protocol)
    np.savez_compressed(out/'declared-inputs.npz',parameters=normalized,initial=initial,scales=np.array(scales),omega=np.array(omegas),
        equilibria=np.array(equilibria),tuple_id=np.repeat(np.arange(len(physical)),2),start_id=np.tile([0,1],len(physical)))
    r=gpu.scan(normalized,initial);arrays={k:v for k,v in r.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(out/'complete-or-budgeted.npz',**arrays);write(out/'settings.json',{k:v for k,v in r.items() if k not in arrays})
    good=r['status']==0
    retained=np.flatnonzero(good&(r['lle_windows'].min(1)>.005));weak=np.flatnonzero(good&(r['lle_windows'].min(1)>.001))
    write(out/'analysis.json',dict(lanes=len(initial),statuses={str(int(k)):int(v) for k,v in zip(*np.unique(r['status'],return_counts=True))},
        retained_ids=retained.tolist(),weak_ids=weak.tolist(),wall_seconds=r['wall_seconds']))
    print('retained',retained,'weak',weak,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--validate-only',action='store_true')
    p.add_argument('--validation',type=Path);main(p.parse_args())

"""Declared exploratory kinetic survey after the restricted ADSI pilot."""
import argparse
import hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import qmc
from search_rho_adsi import GPU,write
from rho_adsi_kinetics import field_jac,independent_field,equilibria,parameterize
from rho_adsi import field_jac as restricted_field


def controls(gpu,out):
    rows=np.array([[1,1.05,1,1,1,1,1,1],[.3,10,.2,4,1,1,1,1],
                   [3,.1,5,.1,1,1,1,1],[1,2,.3,.7,.2,2,.5,3]])
    parameters=parameterize(rows);records=[]
    for i,k in enumerate(parameters):
        x=np.array([.5,3.,1.]);v=np.array([1.,2.,3.])/np.sqrt(14.)
        f,j=field_jac(x,k)
        j2=np.column_stack([independent_field(x+1e-25j*np.eye(3)[a],k).imag/1e-25 for a in range(3)])
        algebra=float(np.max(abs(j-j2)/(1+abs(j2))))
        assert algebra<1e-12
        if i<3:
            r=rows[i];p=[r[0],r[0]*r[1],r[2],r[3]]
            oldf,oldj=restricted_field(x,p)
            assert np.max(abs(f-oldf)/(1+abs(oldf)))<1e-12
            assert np.max(abs(j-oldj)/(1+abs(oldj)))<1e-12
        def rhs(t,z):
            f,j=field_jac(z[:3],k);return np.r_[f,j@z[3:]]
        sol=solve_ivp(rhs,(0,2),np.r_[x,v],method='Radau',rtol=1e-11,atol=1e-13,max_step=.002)
        assert sol.success
        result=gpu.scan(k[None],x[None],0,2,rtol=1e-11,maxstep=.005)
        norm=np.linalg.norm(sol.y[3:,-1])
        errors=dict(algebra=algebra,state=float(np.max(abs(result['states'][0,:3]-sol.y[:3,-1]))),
                    tangent=float(np.max(abs(result['states'][0,3:6]-sol.y[3:,-1]/norm))),
                    rate=float(abs(result['rate'][0]-np.log(norm)/2)))
        assert result['status'][0]==0 and max(errors.values())<1e-7
        records.append(dict(parameters=k.tolist(),errors=errors))
    write(out/'independent-controls.json',records)


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    low=np.array([.1,.5,.005,.005,.01,.2,.05,.05]);high=np.array([10,8,20,20,100,5,20,20])
    seed=2026103147;u=qmc.Sobol(8,scramble=True,seed=seed).random_base2(args.power)
    design=np.exp(np.log(low)+u*np.log(high/low));allk=parameterize(design)
    keep=[];eq=[];rejected_asymptotic=0;residual=0.
    for i,k in enumerate(allk):
        high_activation=k[1]/k[2]*(k[5]/k[6])
        high_inhibition=k[3]+k[4]*(k[7]+k[8]/k[9])/k[10]
        if high_activation>=high_inhibition:
            rejected_asymptotic+=1;continue
        for root in equilibria(k):
            if root['state'][0]<20 and max(root['eigenvalues'].real)>.002:
                keep.append(i);eq.append(root['state']);residual=max(residual,root['residual']);break
    assert len(keep)>0 and residual<1e-7
    params=np.repeat(allk[keep],2,axis=0)
    starts=np.empty((len(params),3));starts[::2]=eq;starts[::2,0]*=1.01
    starts[1::2]=np.column_stack([np.full(len(keep),.5),allk[keep,5]/allk[keep,6],np.ones(len(keep))])
    write(out/'protocol.json',dict(source='Chomchai et al.2025 Eq7/TableS2; same reaction laws, varied kinetic constants.',
        seed=seed,proposals=len(design),selected_tuples=len(keep),lanes=len(params),lower=low.tolist(),upper=high.tolist(),
        design_parameters=['GEF','GAP/GEF','inhibitor_speed','exchange_speed','k0_factor','inactive_pool_factor','k2_factor','k9_factor'],
        parameter_storage='First11 entries: k0..k10 with GEF/GAP folded into k1/k4; last5 padding zeros.',
        selected='high_activation<high_inhibition; positive equilibrium with RT<20, maxRe(eigenvalue)>.002/tau; first such numerical root.',
        high_activation='(k1/k2)*(k5/k6)',high_inhibition='k3+k4*(k7+k8/k9)/k10',
        selection_scope='Heuristic only; not a global boundedness proof or exclusion of hidden attractors.',
        rejected_asymptotic_heuristic=rejected_asymptotic,equilibrium_residual_max=residual,
        initial='1% active perturbation of selected equilibrium and (.5,k5/k6,1).',
        warm=args.warm,observation=args.duration,rtol=1e-9,atol=1e-12,maxstep=.05,
        time='tau=.025*t; physical t in seconds',retention='All4 leading windows>.005/tau; separately declared exploratory>.001.',
        limitation='Broad exploratory source-parameter family; no experimental calibration or novelty claim.'))
    np.savez_compressed(out/'declared-inputs.npz',design=design,parameters=params,initial=starts,
                        selected_proposal_ids=keep,selected_equilibria=eq)
    for name in ['search_rho_adsi_kinetics.py','rho_adsi_kinetics.py','rho_adsi_kinetics.cl','search_rho_adsi.py','rho_adsi.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('declared',len(params),'lanes',flush=True)
    gpu=GPU(kernel='rho_adsi_kinetics.cl',field=field_jac);controls(gpu,out)
    result=gpu.scan(params,starts,args.warm,args.duration,wall=args.wall)
    np.savez_compressed(out/'complete-or-budgeted.npz',**result)
    primary=(result['status']==0)&(result['windows'].min(axis=1)>.005)
    weak=(result['status']==0)&(result['windows'].min(axis=1)>.001)
    write(out/'analysis.json',dict(total=len(params),statuses={str(i):int((result['status']==i).sum()) for i in range(6)},
        primary_ids=np.flatnonzero(primary).tolist(),primary_rates=result['rate'][primary].tolist(),
        exploratory_weak_ids=np.flatnonzero(weak).tolist(),exploratory_weak_rates=result['rate'][weak].tolist(),
        wall_seconds=result['wall_seconds'],initial_field_error=result['initial_field_error'],
        sha256=hashlib.sha256((out/'complete-or-budgeted.npz').read_bytes()).hexdigest()))
    print('completed',int((result['status']==0).sum()),'primary',int(primary.sum()),'weak',int(weak.sum()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--power',type=int,default=15);p.add_argument('--warm',type=float,default=500.)
    p.add_argument('--duration',type=float,default=2000.);p.add_argument('--wall',type=float,default=1200.)
    main(p.parse_args())

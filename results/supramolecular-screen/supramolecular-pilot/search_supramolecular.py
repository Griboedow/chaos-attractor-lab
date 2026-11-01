"""Declared search in the published approximate polymer moment closure."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import qmc
from search_rho_adsi import GPU,write
from supramolecular import REFERENCE,field_jac,independent_field,equilibrium


def controls(gpu,out):
    records=[]
    for p in [REFERENCE,np.array([1.,200.,8.,4.6]),np.array([.03,5000.,3.,4.6])]:
        x=equilibrium(p)['state']*np.array([1.001,.999,1.001]);v=np.array([1.,2.,3.])/np.sqrt(14.)
        f,j=field_jac(x,p);jj=np.column_stack([independent_field(x+1e-25j*np.eye(3)[i],p).imag/1e-25 for i in range(3)])
        algebra=float(np.max(abs(j-jj)/(1+abs(jj))));assert algebra<1e-10
        def augmented(t,z):
            f,j=field_jac(z[:3],p);return np.r_[f,j@z[3:]]
        cpu=solve_ivp(augmented,(0,2),np.r_[x,v],method='Radau',rtol=1e-11,atol=1e-13,max_step=.002)
        assert cpu.success
        result=gpu.scan(p[None],x[None],0,2,rtol=1e-11,maxstep=.005)
        norm=np.linalg.norm(cpu.y[3:,-1]);errors=dict(algebra=algebra,state=float(np.max(abs(result['states'][0,:3]-cpu.y[:3,-1]))),
            tangent=float(np.max(abs(result['states'][0,3:6]-cpu.y[3:,-1]/norm))),rate=float(abs(result['rate'][0]-np.log(norm)/2)))
        assert result['status'][0]==0 and max(errors.values())<1e-7
        records.append(dict(parameters=p.tolist(),errors=errors))
    write(out/'independent-controls.json',records)


def main(args):
    control=json.loads((args.source_control/'analysis.json').read_text());assert control['independent_four_state_solver_scaled_error']<1e-6
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    low=np.array([.005,10.,1.05]);high=np.array([20.,1e4,100.]);seed=2026103149
    u=qmc.Sobol(3,scramble=True,seed=seed).random_base2(args.power)
    design=np.column_stack([np.exp(np.log(low)+u*np.log(high/low)),np.full(len(u),4.6)])
    selected=[REFERENCE];selected_ids=[-1];eqs=[equilibrium(REFERENCE)['state']]
    for i,p in enumerate(design):
        eq=equilibrium(p)
        if max(eq['eigenvalues'].real)>.001 and eq['length']>10:
            selected.append(p);selected_ids.append(i);eqs.append(eq['state'])
    params=np.repeat(selected,2,axis=0);starts=np.empty((len(params),3))
    for i,(p,eq) in enumerate(zip(selected,eqs)):
        # Start from physical mass allocations, not arbitrary state-coordinate perturbations.
        mass=p[2]-eq[0]-eq[1]
        for j,(da,dn) in enumerate([(.01,.01),(-.01,-.01)]):
            x=eq[0]+da*min(eq[0],mass);n=eq[2]*(1+dn);y=eq[1]
            starts[2*i+j]=[x,y,n];assert p[1]*(p[2]-x-y)>2*n
    write(out/'protocol.json',dict(source='Sharko2022 simplified moment closure, thesis Eqs5.6/5.12; Singh2026 Eqs1-4.',
        source_control_sha256=hashlib.sha256((args.source_control/'analysis.json').read_bytes()).hexdigest(),
        seed=seed,proposals=len(design),selected_tuples=len(selected),lanes=len(params),lower=low.tolist(),upper=high.tolist(),
        parameters=['alpha=ka/kd3','beta=kd2/kd3','T=ctot/kd2','C=4.6 fixed'],state='x=a1/kd2,y=d/kd2,z=m0/kd3; k+=1.',
        selection='Source plus positive equilibrium with maxRe>.001 and mean length>10. Heuristic, not a global exclusion.',
        initial='Two opposite 1% active-monomer/mass exchanges and polymer-number perturbations around the equilibrium; declared exact arrays.',
        warm=100.,observation=800.,rtol=1e-9,atol=1e-12,maxstep=.02,time='tau=kd3*t, original source t in units1/(k+K).',
        retention='All4 leading windows>.005; separately exploratory>.001. Minimum observed polymer length>10 required for long-polymer-closure interpretation.',
        error_norm='Componentwise state tolerance; tangent tolerance scaled by full tangent norm.',
        limitations='Approximate long-polymer moment closure, not full size network. No new-model, calibration, chaos or historical-priority claim.'))
    np.savez_compressed(out/'declared-inputs.npz',design=design,selected_proposal_ids=selected_ids,selected_equilibria=eqs,parameters=params,initial=starts)
    for name in ['search_supramolecular.py','supramolecular.py','supramolecular.cl','search_rho_adsi.py','rho_adsi.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('declared',len(params),'lanes',flush=True);gpu=GPU(kernel='supramolecular.cl',field=field_jac);controls(gpu,out)
    result=gpu.scan(params,starts,100,800,wall=1200,maxstep=.02)
    np.savez_compressed(out/'complete-or-budgeted.npz',**result)
    primary=(result['status']==0)&(result['windows'].min(axis=1)>.005);valid=primary&(result['stats'][:,7]>10)
    weak=(result['status']==0)&(result['windows'].min(axis=1)>.001)
    write(out/'analysis.json',dict(total=len(params),statuses={str(i):int((result['status']==i).sum()) for i in range(6)},
        primary_ids=np.flatnonzero(primary).tolist(),long_polymer_primary_ids=np.flatnonzero(valid).tolist(),
        exploratory_ids=np.flatnonzero(weak).tolist(),wall_seconds=result['wall_seconds'],initial_field_error=result['initial_field_error'],
        sha256=hashlib.sha256((out/'complete-or-budgeted.npz').read_bytes()).hexdigest()))
    print('complete',int((result['status']==0).sum()),'primary',int(primary.sum()),'long_polymer_primary',int(valid.sum()),'weak',int(weak.sum()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-control',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--power',type=int,default=13);main(p.parse_args())

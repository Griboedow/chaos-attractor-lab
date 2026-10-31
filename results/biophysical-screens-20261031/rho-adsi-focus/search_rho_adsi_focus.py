"""ADSI continuation near unstable positive equilibria; declared adaptive search."""
import argparse
import hashlib
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from rho_adsi import equilibria
from search_rho_adsi import GPU,write,controls


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    seed=2026103144;u=qmc.Sobol(4,scramble=True,seed=seed).random_base2(args.power)
    low=np.array([.1,.8,.02,.02]);high=np.array([10.,5.,10.,30.])
    raw=np.exp(np.log(low)+u*np.log(high/low));p=raw.copy();p[:,1]*=p[:,0]
    selected=[];eq=[];proposal_ids=[];residual=0.
    for i,row in enumerate(p):
        for root in equilibria(row):
            if root['state'][0]<20 and max(root['eigenvalues'].real)>.002:
                selected.append(row);eq.append(root['state']);proposal_ids.append(i)
                residual=max(residual,root['residual']);break
    params=np.repeat(selected,2,axis=0)
    starts=np.tile([.5,3.,1.],(len(params),1));starts[::2]=np.array(eq);starts[::2,0]*=1.01
    assert residual<1e-8
    write(out/'protocol.json',dict(seed=seed,proposals=len(p),selected_tuples=len(selected),lanes=len(params),
        design='Scrambled Sobol log parameters GEF,GAP/GEF,inhibitor_speed,exchange_speed',
        low=low.tolist(),high=high.tolist(),selection='At least one numerical positive equilibrium with RT<20 and max Re(lambda)>.002 per tau; first such root, perturbed RT by1%, plus fixed (.5,3,1) start.',
        exclusion_scope='A search heuristic, not a proof that excluded tuples lack hidden attractors or that retained trajectories are bounded.',
        warm=args.warm,observation=args.duration,method='Adaptive FP64 DP5(4) state and analytic tangent; rtol1e-9,atol1e-12,maxstep.05tau',
        time='tau=.025*t, t in seconds',equilibrium_residual_max=residual,
        rationale='Follow-up after broad pilot; concentrate on local instability and avoid most large-concentration drift cases.',
        warning='No extra reactions, clipping, forcing or model fitting. Not experimental calibration.'))
    np.savez_compressed(out/'declared-inputs.npz',parameters=params,initial=starts,proposal_parameters=p,
                        selected_proposal_ids=proposal_ids,selected_equilibria=eq)
    for name in ['search_rho_adsi_focus.py','search_rho_adsi.py','rho_adsi.py','rho_adsi.cl']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('declared',len(params),'lanes',flush=True)
    gpu=GPU();controls(gpu,out)
    result=gpu.scan(params,starts,args.warm,args.duration,wall=args.wall)
    np.savez_compressed(out/'complete-or-budgeted.npz',**result)
    primary=(result['status']==0)&(result['windows'].min(axis=1)>.005)
    weak=(result['status']==0)&(result['windows'].min(axis=1)>.001)
    write(out/'analysis.json',dict(total=len(params),statuses={str(i):int((result['status']==i).sum()) for i in range(6)},
        primary_ids=np.flatnonzero(primary).tolist(),primary_rates=result['rate'][primary].tolist(),
        exploratory_weak_ids=np.flatnonzero(weak).tolist(),exploratory_weak_rates=result['rate'][weak].tolist(),
        primary='all four leading windows>.005/tau',exploratory='all four leading windows>.001/tau, explicitly lower threshold',
        wall_seconds=result['wall_seconds'],initial_field_error=result['initial_field_error'],
        sha256=hashlib.sha256((out/'complete-or-budgeted.npz').read_bytes()).hexdigest()))
    print('completed',int((result['status']==0).sum()),'primary',int(primary.sum()),'weak',int(weak.sum()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--power',type=int,default=14);p.add_argument('--warm',type=float,default=1000.)
    p.add_argument('--duration',type=float,default=4000.);p.add_argument('--wall',type=float,default=1200.)
    main(p.parse_args())

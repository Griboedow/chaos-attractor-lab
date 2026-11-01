"""Test all row permutations against phase shifts of a numerical motor cycle."""
import argparse
import hashlib
from itertools import permutations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from audit_motor_permutation_geometry import full
from motor_rows import rhs


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    data=np.load(args.source);seed=data['seed'];p=data['parameters'];period=float(data['period'])
    sol=solve_ivp(lambda t,q:rhs(t,q,p,4),(0,period),seed,method='DOP853',rtol=2e-12,atol=2e-14,max_step=.005,dense_output=True)
    assert sol.success
    weights=np.r_[np.ones(8),np.full(4,p[0]/p[1])]/np.sqrt(1+p[2]**2)
    grid=np.linspace(0,period,4096,endpoint=False);orbit=full(sol.sol(grid).T,4);dt=period/len(grid)
    records=[]
    # Optimize one seed, then verify the same shift at 32 independent phases.
    validation=np.linspace(0,period,32,endpoint=False)
    original=full(sol.sol(validation).T,4)
    for perm in permutations(range(4)):
        ix=np.r_[perm,np.array(perm)+4,np.array(perm)+8];target=full(seed,4)[ix]
        error=np.linalg.norm((orbit-target)*weights,axis=1);k=int(np.argmin(error))
        def objective(t):return float(np.linalg.norm((full(sol.sol(t%period),4)-target)*weights)**2)
        optimum=minimize_scalar(objective,bounds=(grid[k]-2*dt,grid[k]+2*dt),method='bounded',options={'xatol':1e-13})
        shift=float(optimum.x%period);shifted=full(sol.sol((validation+shift)%period).T,4)
        residual=(shifted-original[:,ix])*weights
        records.append(dict(permutation=list(perm),phase_shift=shift,phase_fraction=shift/period,
            seed_distance=float(np.sqrt(optimum.fun)),maximum_32phase_error=float(np.max(abs(residual))),
            close_at_1e_8=bool(np.max(abs(residual))<1e-8)))
    report=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),parameters=p.tolist(),
        period=period,closure=float(np.max(abs(sol.y[:,-1]-seed))),permutations=records,
        method='All24 exact row permutations;4096-phase coarse search, scalar bounded refinement,32-phase full12D verification.',
        limitation='Finite floating-point spatiotemporal symmetry test, not a proof of group isotropy or attracting-set identity.')
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    for name in ['audit_motor_cycle_symmetry.py','audit_motor_permutation_geometry.py','motor_rows.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print([(r['permutation'],r['phase_fraction'],r['maximum_32phase_error']) for r in records if r['close_at_1e_8']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

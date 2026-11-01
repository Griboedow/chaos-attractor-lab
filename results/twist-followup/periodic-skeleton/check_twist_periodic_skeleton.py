"""Numerical fixed points and three-cycles of the archived twist geometries.

This is a non-exhaustive root search, not a proof of chaos or of its mechanism.
Every root start, solver status and residual is retained.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.spatial import cKDTree
from scipy.stats import qmc
from twist_map import step,step_tangent


def iterate(z,p,period):
    q=np.eye(3)
    for _ in range(period):z,q=step_tangent(z,q,p)
    return z,q


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=args.source;long=np.load(source/'long-runs.npz');basins=np.load(source/'basins.npz')
    selected=[1709,3254,307];seed=2026110116
    cube=2*qmc.Sobol(3,scramble=True,seed=seed).random_base2(11)-1
    sphere=cube[np.linalg.norm(cube,axis=1)<=1]
    protocol=dict(sources={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in ['long-runs.npz','basins.npz','candidates.json']},
        candidates=selected,sobol_seed=seed,fixed_point_starts='All unit-ball points in a 2048-row Sobol cube, plus 128 uniformly indexed reference-tail states.',
        period_three_starts='256 reference-tail states with smallest three-step return distance, plus one last state from every archived fresh-start run.',
        root='SciPy hybr with analytic derivative,xtol1e-11. Accept residual<1e-10 and radius<=1+1e-8 regardless of solver status; retain every status.',
        deduplication='Orbit points agree within1e-7, allowing cyclic shifts. Check period1 first, then3.',
        scope='Non-exhaustive numerical periodic-orbit skeleton. No global classification, homoclinic connection, attractor membership or novelty proof.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_twist_periodic_skeleton.py','twist_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    for ident in selected:
        index=list(long['ids']).index(ident);p=long['parameters'][index];tail=long['tail'][index];tree=cKDTree(tail)
        fixed_starts=np.r_[sphere,tail[np.linspace(0,len(tail)-1,128).astype(int)]]
        return_order=np.argsort(np.linalg.norm(tail[3:]-tail[:-3],axis=1))[:256]
        three_starts=np.r_[tail[return_order],basins['tail'][basins['ids']==ident,-1]]
        attempts=[];orbits=[]
        for period,starts in [(1,fixed_starts),(3,three_starts)]:
            for initial in starts:
                sol=root(lambda z:iterate(z,p,period)[0]-z,initial,
                         jac=lambda z:iterate(z,p,period)[1]-np.eye(3),method='hybr',tol=1e-11)
                z=sol.x;residual=float(np.linalg.norm(iterate(z,p,period)[0]-z));radius=float(np.linalg.norm(z))
                accepted=bool(np.isfinite(z).all() and residual<1e-10 and radius<=1+1e-8)
                attempts.append(dict(period=period,initial=initial.tolist(),status=int(sol.status),nfev=int(sol.nfev),
                    state=z.tolist(),residual=residual,radius=radius,accepted=accepted))
                if not accepted:continue
                least=1 if np.linalg.norm(step(z,p)-z)<1e-9 else period
                orbit=[z.copy()]
                for _ in range(least-1):orbit.append(step(orbit[-1],p))
                orbit=np.array(orbit)
                if any(a['least_period']==least and np.min(np.linalg.norm(np.array(a['states'])-z,axis=1))<1e-7 for a in orbits):continue
                end,j=iterate(z,p,least);eig=np.linalg.eigvals(j)
                distances=tree.query(orbit)[0]
                orbits.append(dict(least_period=least,states=orbit.tolist(),residual=float(np.linalg.norm(end-z)),
                    multipliers=[[float(a.real),float(a.imag)] for a in eig],moduli=abs(eig).tolist(),
                    unstable_directions=int(np.count_nonzero(abs(eig)>1+1e-7)),
                    volume_identity_error=float(abs(np.linalg.det(j)-p[0]**(3*least))),
                    nearest_reference_distance=distances.tolist()))
            print(ident,'period',period,'orbits',len(orbits),flush=True)
        row=dict(id=ident,parameters=p.tolist(),attempts=len(attempts),accepted_attempts=sum(a['accepted'] for a in attempts),
                 periodic_orbits=orbits,reference_samples=len(tail))
        (out/f'attempts-{ident}.json').write_text(json.dumps(attempts,indent=2,allow_nan=False)+'\n')
        records.append(row);(out/'analysis.json').write_text(json.dumps(dict(runs=records),indent=2,allow_nan=False)+'\n')
        print('summary',ident,[(a['least_period'],a['unstable_directions'],min(a['nearest_reference_distance'])) for a in orbits],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

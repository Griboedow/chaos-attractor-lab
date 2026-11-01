"""Finite diagnostics for posterior switching, folds and singular preimages.

These calculations propose geometry and obstructions. They do not certify
trapping, invariant membership, physical measures, or historical novelty.
"""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from numba import njit
from scipy.optimize import root
from scipy.spatial import cKDTree
from scipy.stats import qmc
from search_product_alphabet import step


@njit(cache=True)
def diagnostics(states, beta, tau, levels):
    count = len(states)
    data = np.empty((count, 8))
    signs = np.empty(count, np.int64)
    vector = np.array([1., .371, -.217]); vector /= np.linalg.norm(vector)
    for k in range(count):
        x = states[k]
        _, J, r, _, entropy, flag = step(x, beta, tau, levels, 2)
        assert not flag
        v = J @ vector; length = np.linalg.norm(v); vector = v / length
        d = np.empty(3); label = 0
        for j in range(3):
            lo, hi = levels[j, 0], levels[j, 1]
            a, c = (hi-lo)/2, (hi+lo)/2
            h = 2*a*(x[j]-c)/tau
            e = np.exp(-2*abs(h))
            d[j] = 8*a*a/tau*e/(1+e)**2
            if x[j] > c: label += 1 << j
        signs[k] = label
        data[k, 0] = r
        data[k, 1] = np.linalg.det(J)
        data[k, 2] = np.log(length)
        data[k, 3:6] = d
        data[k, 6] = np.sqrt(np.sum(J*J))
        data[k, 7] = entropy
    return data, signs


@njit(cache=True)
def images(states, beta, tau, levels):
    out = np.empty_like(states)
    for k in range(len(states)):
        out[k], _, _, _, _, flag = step(states[k], beta, tau, levels, 2)
        assert not flag
    return out


def main(args):
    args.out.mkdir(parents=True, exist_ok=False); begin = time.time()
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    names = [Path(__file__).name, 'search_product_alphabet.py', 'follow_spherical_constellations.py',
             'search_soft_sphere.py', 'search_soft_constellations.py', 'search_raar_constellations.py',
             'search_spherical_constellations.py']
    for name in names: (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:
        states = f['states']; beta, tau = f['params']; levels = f['levels'][:, :2].copy()
    d, signs = diagnostics(states, beta, tau, levels)
    transitions = np.bincount(signs[:-1]*8+signs[1:], minlength=64).reshape(8,8)
    active = (d[:,3:6]>1).sum(1)
    groups = []
    for n in range(4):
        mask = (active == n); mask[:1000] = False
        groups.append(dict(active_coordinates=n,count=int(mask.sum()),
                           mean_leading_log=float(d[mask,2].mean()) if mask.any() else None,
                           leading_log_sum=float(d[mask,2].sum()),
                           negative_determinant=int((d[mask,1]<0).sum())))
    c = levels.mean(1)
    def singular_equation(x):
        a = np.diff(levels, axis=1)[:,0]/2
        m = c + a*np.tanh(2*a*(x-c)/tau)
        return 2*m-x
    singular = root(singular_equation, c, tol=1e-12).x
    assert np.max(abs(singular_equation(singular)))<1e-12
    tree = cKDTree(states)
    ys = images(states[::20], beta,tau,levels)
    nearest = np.argsort(np.linalg.norm(ys-singular,axis=1))[:64]
    starts = np.vstack([states[::20][nearest],1.6*qmc.Sobol(3,scramble=True,seed=2026110156).random_base2(9)-.8])
    preimages, attempts = [], []
    for initial in starts:
        fit = root(lambda x: step(x,beta,tau,levels,2)[0]-singular, initial,
                   jac=lambda x: step(x,beta,tau,levels,2)[1], tol=1e-11)
        err = float(np.max(abs(step(fit.x,beta,tau,levels,2)[0]-singular)))
        found = None
        if err<1e-10 and np.isfinite(fit.x).all():
            found = next((k for k,p in enumerate(preimages) if np.linalg.norm(fit.x-p['state'])<1e-7), None)
            if found is None:
                found = len(preimages)
                preimages.append(dict(state=fit.x.tolist(), residual=err,
                                      cloud_distance=float(tree.query(fit.x)[0]),
                                      in_coordinate_box=bool(np.all(fit.x>=states.min(0)) and np.all(fit.x<=states.max(0))),
                                      determinant=float(np.linalg.det(step(fit.x,beta,tau,levels,2)[1]))))
        attempts.append(dict(start=initial.tolist(), success=bool(fit.success), residual=err if np.isfinite(err) else None, preimage=found))
    cycles = json.loads(args.cycles.read_text(encoding='utf-8'))['cycles']
    cycle_summary = [dict(id=p['id'],period=p['period'],unstable_directions=p['unstable_directions'],
                         maximum_cloud_distance=float(tree.query(np.array(p['orbit']))[0].max()), multipliers=p['multipliers']) for p in cycles]
    np.savez_compressed(args.out/'diagnostics.npz', values=d, signs=signs, params=[beta,tau], levels=levels)
    report = dict(source_sha256=sha(args.source), source=str(args.source), cycles_source_sha256=sha(args.cycles),
                  states=len(states), params=[float(beta),float(tau)], levels=levels.tolist(),
                  coordinate_min=states.min(0).tolist(), coordinate_max=states.max(0).tolist(),
                  sign_counts=np.bincount(signs,minlength=8).tolist(), transitions=transitions.tolist(),
                  conditional_stretch=groups, leading_mean=float(d[1000:,2].mean()),
                  radius_quantiles=np.quantile(d[:,0],[0,.0001,.001,.01,.5,.99,1]).tolist(),
                  determinant_quantiles=np.quantile(d[:,1],[0,.001,.01,.5,.99,.999,1]).tolist(),
                  minimum_cloud_distance_to_central_singularity=float(tree.query(singular)[0]),
                  central_singularity_numerical=singular.tolist(), preimages=preimages,
                  preimage_attempts=attempts,cycles=cycle_summary,seconds=time.time()-begin,
                  scope=__doc__)
    (args.out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report[k] for k in ['coordinate_min','coordinate_max','sign_counts','conditional_stretch','leading_mean','radius_quantiles','minimum_cloud_distance_to_central_singularity','preimages','seconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','cycles','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())

"""Nonuniform grid-closure proposals. Sampled remainders are NOT proofs."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.optimize import brentq
from search_product_alphabet import step


@njit(cache=True)
def occupied(mask):
    return np.flatnonzero(mask)


@njit(cache=True)
def advance(frontier, mask, edges, shape, beta, tau, levels, order, curvature, slack, budget):
    ny, nz = shape[1], shape[2]
    added = [np.int64(0) for _ in range(0)]
    minr, maxrem = 1e100, 0.
    for pos, key in enumerate(frontier):
        ix = np.array([key//(ny*nz), (key//nz)%ny, key%nz])
        mid, rad = np.zeros(3), np.zeros(3)
        for i in range(3):
            mid[i] = (edges[i, ix[i]]+edges[i, ix[i]+1])/2
            rad[i] = (edges[i, ix[i]+1]-edges[i, ix[i]])/2
        center, J, r, _, _, flag = step(mid, beta, tau, levels, order)
        minr = min(minr, r)
        if flag: return np.array(added), 1, pos, key, minr, maxrem
        rem = np.zeros(3)
        for c in range(8):
            delta = rad*np.array([2*(c&1)-1, 2*((c>>1)&1)-1, 2*((c>>2)&1)-1])
            y, _, r, _, _, flag = step(mid+delta, beta, tau, levels, order)
            minr = min(minr, r)
            if flag: return np.array(added), 1, pos, key, minr, maxrem
            rem = np.maximum(rem, np.abs(y-center-J@delta))
        rem = curvature*rem+slack
        maxrem = max(maxrem, np.max(rem))
        radius = np.abs(J)@rad+rem
        lower, upper = np.zeros(3, np.int64), np.zeros(3, np.int64)
        for i in range(3):
            lower[i] = np.searchsorted(edges[i, :shape[i]+1], center[i]-radius[i], side='right')-1
            upper[i] = np.searchsorted(edges[i, :shape[i]+1], center[i]+radius[i], side='right')-1
        if np.any(lower<0) or np.any(upper>=shape):
            return np.array(added), 2, pos, key, minr, maxrem
        directions = np.zeros((15, 3)); directions[:3] = np.eye(3); k = 3
        for i in range(3):
            for j in range(i+1, 3): directions[k] = np.cross(J[:, i], J[:, j]); k += 1
        for i in range(3):
            for j in range(3): directions[k] = np.cross(J[:, i], np.eye(3)[j]); k += 1
        support = np.zeros(15)
        for k in range(15):
            length = np.linalg.norm(directions[k])
            if length>1e-15: directions[k] /= length
            support[k] = np.sum(np.abs(directions[k]@J)*rad)+np.sum(np.abs(directions[k])*rem)
        for a in range(lower[0], upper[0]+1):
            for b in range(lower[1], upper[1]+1):
                for c in range(lower[2], upper[2]+1):
                    target = (a*ny+b)*nz+c
                    if mask[target]: continue
                    ti = np.array([a, b, c]); displacement, tr = np.zeros(3), np.zeros(3)
                    for i in range(3):
                        displacement[i] = (edges[i, ti[i]]+edges[i, ti[i]+1])/2-center[i]
                        tr[i] = (edges[i, ti[i]+1]-edges[i, ti[i]])/2
                    intersects = True
                    for k in range(15):
                        if np.abs(np.dot(directions[k], displacement)) > support[k]+np.sum(np.abs(directions[k])*tr):
                            intersects = False; break
                    if intersects: mask[target] = True; added.append(target)
        if len(added)>budget: return np.array(added), 3, pos, key, minr, maxrem
    return np.array(added), 0, len(frontier), -1, minr, maxrem


def main(args):
    args.out.mkdir(parents=True, exist_ok=False); start = time.time()
    for name in [Path(__file__).name, 'search_product_alphabet.py', 'search_soft_sphere.py', 'search_soft_constellations.py', 'search_raar_constellations.py', 'follow_spherical_constellations.py', 'search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:
        z, params, levels, order = f['states'], f['params'], f['levels'], int(f['order'])
    assert order == 2
    beta, tau = params; lists = []
    for i in range(3):
        c, a = np.mean(levels[i, :2]), np.diff(levels[i, :2])[0]/2
        def coordinate(x): return x+args.alpha*(c+a*np.tanh(2*a*(x-c)/tau))
        low, high = -args.halfwidth, args.halfwidth
        ys = np.linspace(coordinate(low), coordinate(high), int(np.ceil((coordinate(high)-coordinate(low))/args.step))+1)
        lists.append(np.array([low]+[brentq(lambda x: coordinate(x)-y, low, high, xtol=1e-15) for y in ys[1:-1]]+[high]))
    shape = np.array([len(e)-1 for e in lists], dtype=np.int64)
    edges = np.full((3, max(shape)+1), np.inf)
    for i in range(3): edges[i, :len(lists[i])] = lists[i]
    assert np.prod(shape)<args.max_cells
    ijk = np.column_stack([np.searchsorted(lists[i], z[:, i], side='right')-1 for i in range(3)])
    assert np.all(ijk>=0) and np.all(ijk<shape)
    frontier = np.unique((ijk[:, 0]*shape[1]+ijk[:, 1])*shape[2]+ijk[:, 2])
    mask = np.zeros(np.prod(shape), dtype=np.bool_); mask[frontier] = True
    initial = len(frontier); history = []; flag = 6; key = -1
    print(json.dumps(dict(shape=shape.tolist(), initial=initial, min_width=[float(np.diff(e).min()) for e in lists], max_width=[float(np.diff(e).max()) for e in lists])), flush=True)
    for generation in range(args.generations):
        new, flag, pos, key, minr, maxrem = advance(frontier, mask, edges, shape, beta, tau, levels, order, args.curvature, args.slack, args.budget)
        row = dict(generation=generation, sources=len(frontier), processed=int(pos), added=len(new), status=int(flag), minimum_sampled_radius=float(minr), maximum_sampled_remainder=float(maxrem), seconds=time.time()-start)
        history.append(row); print(json.dumps(row), flush=True)
        (args.out/'progress.json').write_text(json.dumps(history, indent=2)+'\n', encoding='utf-8', newline='\n')
        if flag: break
        frontier = new
        if not len(frontier): break
    keys = occupied(mask)
    np.savez_compressed(args.out/'proposal.npz', edges=edges, shape=shape, occupied=keys, pending=frontier, params=params, levels=levels, order=order)
    result = dict(status='numerically_closed' if flag==0 and not len(frontier) else 'not_closed', reason={0:'Sampled closure or generation limit', 1:'Undefined sampled map', 2:'Enclosure leaves domain', 3:'New-cell budget exceeded', 6:'No generations'}[int(flag)], source=str(args.source), source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(), arguments={k:str(v) if isinstance(v, Path) else v for k,v in vars(args).items()}, shape=shape.tolist(), initial_occupied=initial, final_occupied=len(keys), failing_key=int(key), history=history, scope='Proposal only: sampled corner residuals are not rigorous bounds. No trapping, attraction, chaos or novelty certificate.')
    (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['history','arguments']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['source', 'out']: p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--step', type=float, default=.02); p.add_argument('--alpha', type=float, default=1.)
    p.add_argument('--halfwidth', type=float, default=.7); p.add_argument('--curvature', type=float, default=2.)
    p.add_argument('--slack', type=float, default=1e-8); p.add_argument('--budget', type=int, default=3000000)
    p.add_argument('--max-cells', type=int, default=100000000); p.add_argument('--generations', type=int, default=100)
    main(p.parse_args())

"""Reject or propose elementary trapping domains using sampled witnesses only."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit, prange, set_num_threads
from scipy.stats import qmc
from search_product_alphabet import step


@njit(cache=True, parallel=True)
def sample(states, beta, tau, levels, order):
    images = np.empty_like(states)
    radius = np.empty((len(states), 2))
    flags = np.empty(len(states), np.int64)
    for i in prange(len(states)):
        y, _, r, _, _, flag = step(states[i], beta, tau, levels, order)
        _, _, s, _, _, other = step(y, beta, tau, levels, order)
        images[i] = y
        radius[i] = r, s
        flags[i] = flag+other
    return images, radius, flags


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    start = time.time()
    set_num_threads(4)
    names = [Path(__file__).name, 'search_product_alphabet.py', 'search_soft_sphere.py', 'search_soft_constellations.py', 'search_raar_constellations.py', 'follow_spherical_constellations.py', 'search_spherical_constellations.py']
    for name in names:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as data:
        params, levels, order = data['params'], data['levels'], int(data['order'])
    states = 2*args.halfwidth*qmc.Sobol(3, scramble=True, seed=args.seed).random_base2(args.power)-args.halfwidth
    images, radius, flags = sample(states, *params, levels, order)
    assert not np.any(flags)
    trials = []
    for kind in ['box', 'ball']:
        for center in [np.zeros(3), np.array([-.05, .086, .022])]:
            distance = np.max(abs(states-center), axis=1) if kind == 'box' else np.linalg.norm(states-center, axis=1)
            output_distance = np.max(abs(images-center), axis=1) if kind == 'box' else np.linalg.norm(images-center, axis=1)
            for outer in args.outer:
                for cutoff in args.cutoff:
                    selected = np.flatnonzero((distance <= outer) & (radius[:, 0] >= cutoff))
                    bad = selected[(output_distance[selected] >= outer) | (radius[selected, 1] <= cutoff)]
                    row = dict(shape=kind, center=center.tolist(), outer=outer, minimum_input_radius=cutoff, sampled_points=len(selected), failures=len(bad), passed_sample=not len(bad), maximum_image_distance=float(output_distance[selected].max()), minimum_image_radius=float(radius[selected, 1].min()))
                    if len(bad):
                        k = int(bad[0])
                        row['witness'] = dict(index=k, state=states[k].tolist(), image=images[k].tolist(), radii=radius[k].tolist())
                    trials.append(row)
    np.savez_compressed(args.out/'sample.npz', states=states, images=images, radius=radius, flags=flags, params=params, levels=levels, order=order)
    result = dict(status='sampled', source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(), seed=args.seed, halfwidth=args.halfwidth, points=len(states), trials=trials, seconds=time.time()-start, scope='Sampled elementary outer domains with a sphere-normalization exclusion. A passing sample is not a trapping certificate. Violating witnesses are finite numerical counterexamples only.')
    (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(trials=len(trials), passing=sum(r['passed_sample'] for r in trials), seconds=result['seconds'])))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['source', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--power', type=int, default=20)
    p.add_argument('--halfwidth', type=float, default=.8)
    p.add_argument('--seed', type=int, default=2026110148)
    p.add_argument('--outer', type=float, nargs='+', default=[.36, .4, .45, .5, .6, .7])
    p.add_argument('--cutoff', type=float, nargs='+', default=[.1, .5, 1., 1.2, 1.5, 2.])
    main(p.parse_args())

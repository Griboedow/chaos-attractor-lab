"""Seeded screening and fixed-step replays of dissipative Cartesian flows."""
import argparse
import csv
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from field import (integrate, make_case, make_ring, make_rotor, make_planes,make_crossed,
                   ring_tangent, rotor_tangent, plane_tangent, velocity_tangent,crossed_tangent)


def metrics(points):
    if points.std(axis=0).max() < 1e-4:
        return [0.] * 6
    eigen = np.linalg.eigvalsh(np.cov(points.T)).clip(0)
    radius = np.linalg.norm(points - points.mean(axis=0), axis=1)
    span = np.ptp(points, axis=0)
    cells = np.minimum(15, (16 * (points - points.min(axis=0)) / np.maximum(span, .001)).astype(int))
    occupied = len(np.unique(np.ravel_multi_index(cells.T, (16, 16, 16))))
    return [round(float(x), 6) for x in (eigen[0] / eigen[-1], eigen[1] / eigen[-1],
            np.quantile(radius, .05) / max(np.quantile(radius, .95), 1e-9), np.linalg.norm(span),
            occupied, np.max(np.abs(points)))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--seed', type=int, default=2026112123)
    parser.add_argument('--family', choices=('cartesian', 'ring', 'rotor', 'planes', 'crossed'), default='cartesian')
    parser.add_argument('--first', type=int, default=0)
    parser.add_argument('--count', type=int, default=256)
    parser.add_argument('--ids', type=int, nargs='+')
    parser.add_argument('--inputs', type=Path, help='Read already selected parameters and initial states')
    parser.add_argument('--dt', type=float, default=.04)
    parser.add_argument('--burn', type=float, default=160.)
    parser.add_argument('--window', type=float, default=80.)
    parser.add_argument('--blocks', type=int, default=4)
    parser.add_argument('--stride', type=int, default=4)
    parser.add_argument('--initial', type=int, default=0)
    parser.add_argument('--local-radius', type=float, help='Displace the original start by this Euclidean distance')
    parser.add_argument('--workers', type=int, choices=(1, 2), default=2)
    parser.add_argument('--save-all', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    ids = args.ids or list(range(args.first, args.first + args.count))
    burn = round(args.burn / args.dt)
    window = round(args.window / args.dt)
    inputs = json.loads(args.inputs.read_text()) if args.inputs else {}

    def run(case):
        generate, rhs = {'cartesian': (make_case, velocity_tangent), 'ring': (make_ring, ring_tangent),
                         'rotor': (make_rotor, rotor_tangent), 'planes': (make_planes, plane_tangent),
                         'crossed': (make_crossed,crossed_tangent)}[args.family]
        p, initial = generate(args.seed, case)
        if inputs:
            p=np.array(inputs[str(case)]['parameters']);initial=np.array(inputs[str(case)]['initial'])
        if args.initial:
            rng = np.random.default_rng(np.random.SeedSequence([args.seed, case, args.initial]))
            if args.local_radius is None:
                initial = rng.uniform(-6., 6., size=3)
            else:
                shift = rng.normal(size=3)
                initial = initial + args.local_radius * shift / np.linalg.norm(shift)
        rates, points = integrate(p, initial, args.dt, burn, window, args.blocks, args.stride, rhs)
        if not np.isfinite(points).all() or not np.isfinite(rates).all():
            return [case, *([float('nan')] * args.blocks), *([0.] * 6), 0]
        shape = metrics(points)
        keep = min(rates) > .008 and shape[0] > .008
        if args.save_all or keep:
            np.savez_compressed(args.out / f'case-{case}.npz', points=points, parameters=p,
                                initial=initial, rates=rates, step_size=args.dt, sample_stride=args.stride)
        return [case, *np.round(rates, 7), *shape, int(keep)]

    (args.out / 'run.json').write_text(json.dumps({**vars(args), 'out': str(args.out)}, indent=2,default=str)+'\n')
    first = run(ids[0])
    with (args.out / 'screen.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['case', *(f'window_{i}' for i in range(args.blocks)), 'variance_ratio_min',
                         'variance_ratio_mid', 'radius_05_over_95', 'span', 'occupied', 'max_abs', 'keep'])
        writer.writerow(first)
        handle.flush()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for i, row in enumerate(pool.map(run, ids[1:]), 2):
                writer.writerow(row)
                handle.flush()
                if i % 32 == 0:
                    print(f'{i}/{len(ids)}', flush=True)


if __name__ == '__main__':
    main()

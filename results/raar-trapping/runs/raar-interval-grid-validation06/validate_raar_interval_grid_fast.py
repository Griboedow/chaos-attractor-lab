"""Exact-rational raster checks and independent high-precision image samples."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

import mpmath as mp
import numpy as np
from check_raar_coverings_ad import Model, exact
from propose_raar_trapping_bitset import mark, marked
import raar_interval_grid as array_grid
import raar_geometric_enclosure as array_geometry
import raar_interval_scalar_grid as scalar_grid
import raar_geometric_scalar as scalar_geometry
import raar_interval_offset_grid_fast as offset_grid
from search_raar_circles import step


def floor_iv(x):
    return [int(mp.floor(mp.mpf(v))) for v in x._mpi_]


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    grid = offset_grid if args.engine == 'offset' else (scalar_grid if args.engine == 'scalar' else array_grid)
    geometry = scalar_geometry if args.engine in ['scalar','offset'] else array_geometry
    offset=np.array([128,128,128] if args.engine=='offset' else [0,0,0],dtype=np.int64)
    advance, coordinate_interval = grid.advance, grid.coordinate_interval
    raster_parameters, target_separated = grid.raster_parameters, grid.target_separated
    remainder = geometry.remainder
    root = Path(__file__).resolve().parents[1]
    names = ['validate_raar_interval_grid_fast.py', 'raar_interval_grid.py', 'raar_interval64.py',
             'raar_geometric_enclosure.py', 'propose_raar_trapping_bitset.py',
             'raar_interval_scalar.py', 'raar_geometric_scalar.py', 'raar_interval_scalar_grid.py','raar_interval_offset_grid_fast.py']
    hashes = {}
    for name in names:
        b = Path(__file__).with_name(name).read_bytes()
        (args.out/name).write_bytes(b)
        hashes[name] = hashlib.sha256(b).hexdigest()
    d = np.load(root/'results/raar-circles/render/case-1051-animation.npz')
    z, p = d['states'], d['params']
    h = .005
    margin = 2.**-35
    origin = np.floor((z.min(0)-.35)/h)*h
    shape = np.ceil((z.max(0)+.35-origin)/h).astype(np.int64)+2*offset
    axis = np.hypot(z[:, 0]-p[3], np.cos(p[2])*(z[:, 1]-p[4])+np.sin(p[2])*(z[:, 2]-p[5]))
    ids = np.unique(np.r_[np.linspace(10000, len(z)-1, 64, dtype=int), np.argsort(axis)[:16]])
    ijk = np.floor((z[ids]-origin)/h).astype(np.int64)+offset
    keys = np.unique((ijk[:, 0]*shape[1]+ijk[:, 1])*shape[2]+ijk[:, 2])
    mask = np.zeros((int(np.prod(shape))+7)//8, np.uint8)
    for key in keys:
        mark(mask, key)
    mp.iv.dps, mp.mp.dps = 100, 120
    trig = np.array([[np.nextafter(float(v.a), -np.inf), np.nextafter(float(v.b), np.inf)]
                     for v in [mp.iv.sin(exact(p[2])), mp.iv.cos(exact(p[2]))]])
    print(json.dumps({'stage': 'compile and rasterize', 'source_cells': len(keys)}), flush=True)
    if args.engine=='offset':
        result=advance(keys,mask,shape,origin,h,p,trig,margin,1000000,4,offset)
        new,flag,processed,failing,minimum,splits,leaves,deepest,bbox_accepts=result;largest=None
    else:
        result = advance(keys, mask, shape, origin, h, p, trig, margin, 1000000, .15, 4)
        new, flag, processed, failing, minimum, splits, leaves, largest = result
    assert flag == 0 and processed == len(keys), (flag, processed, failing)
    rng = np.random.default_rng(202611011046)
    model = Model(p)
    image_checks = rational_separation_checks = 0
    for key in keys:
        index = np.array([key//(shape[1]*shape[2]), (key//shape[2])%shape[1], key%shape[2]])
        box = np.array([[coordinate_interval(origin[j], h, float(index[j]-offset[j]))[0],
                         coordinate_interval(origin[j], h, float(index[j]-offset[j]+1))[1]] for j in range(3)])
        for j in range(3):
            assert Q(float(box[j, 0])) <= Q(float(origin[j]))+Q(h)*int(index[j]-offset[j])
            assert Q(float(box[j, 1])) >= Q(float(origin[j]))+Q(h)*(int(index[j]-offset[j])+1)
        mid = origin+h*(index-offset+.5)
        center, J, _, _ = step(mid, p)
        radii, errors, _, status = remainder(box, mid, center, J, p, trig, margin)
        assert status == 0
        a, b, directions, support = raster_parameters(center,J,radii,errors,origin,h,offset) if args.engine=='offset' else raster_parameters(center, J, radii, errors, origin, h)
        for j in range(3):
            radius = sum(abs(Q(float(J[j, k])))*Q(float(radii[k])) for k in range(3))+Q(float(errors[j]))
            lo = (Q(float(center[j]))-radius-Q(float(origin[j])))/Q(h)
            hi = (Q(float(center[j]))+radius-Q(float(origin[j])))/Q(h)
            assert int(a[j]) <= lo.numerator//lo.denominator+offset[j]
            assert int(b[j]) >= hi.numerator//hi.denominator+offset[j]
        for _ in range(12):
            target = np.array([rng.integers(a[j]-2, b[j]+3) for j in range(3)])
            if (target_separated(target,origin,h,center,directions,support,offset) if args.engine=='offset' else target_separated(target, origin, h, center, directions, support)):
                exact_separated = False
                for normal in directions:
                    n = [Q(float(v)) for v in normal]
                    bound = sum(Q(float(radii[j]))*abs(sum(n[i]*Q(float(J[i, j])) for i in range(3))) for j in range(3))
                    bound += sum(abs(n[i])*(Q(float(errors[i]))+Q(h)/2) for i in range(3))
                    projection = sum(n[j]*(Q(float(origin[j]))+Q(h)*(int(target[j]-offset[j])+Q(1, 2))-Q(float(center[j]))) for j in range(3))
                    exact_separated |= abs(projection) > bound
                assert exact_separated
            rational_separation_checks += 1
        fractions = list(itertools.product([0., 1.], repeat=3))+[tuple(v) for v in rng.random((4, 3))]
        for fraction in fractions:
            x = [exact(origin[j])+exact(h)*(int(index[j]-offset[j])+exact(fraction[j])) for j in range(3)]
            image, _ = model.field(x)
            for sign in itertools.product([-1, 1], repeat=3):
                bounds = [[v+int(offset[j]) for v in floor_iv((image[j]+sign[j]*exact(margin)-exact(origin[j]))/exact(h))] for j in range(3)]
                for target in itertools.product(*[range(lo, hi+1) for lo, hi in bounds]):
                    assert all(0 <= target[j] < shape[j] for j in range(3))
                    target_key = (target[0]*shape[1]+target[1])*shape[2]+target[2]
                    assert marked(mask, target_key), (int(key), target)
                    image_checks += 1
    if args.compare:
        previous = np.load(args.compare/'sample.npz')
        assert np.array_equal(keys, previous['source_keys'])
        assert np.array_equal(np.sort(new), np.sort(previous['new_keys']))
    report = {'engine': args.engine, 'comparison': str(args.compare) if args.compare else None,
              'source_cells': len(keys), 'new_cells': len(new), 'processed': int(processed),
              'subdivision_splits': int(splits), 'accepted_leaves': int(leaves),
              'positive_axis_lower_bound': float(minimum), 'largest_leaf_remainder_over_step': float(largest) if largest is not None else None,
              'exact_rational_separation_checks': rational_separation_checks,
              'independent_high_precision_image_margin_checks': image_checks,
              'source_sha256': hashes, 'seconds': time.time()-started,
              'scope': 'Implementation checks of outward rasterization on selected cells, not a full forward closure or attraction theorem.'}
    np.savez_compressed(args.out/'sample.npz', source_keys=keys, new_keys=new, params=p, trig=trig,
                        origin=origin, step=h, shape=shape, image_margin=margin,grid_offset=offset)
    (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--engine', choices=['array', 'scalar','offset'], default='array')
    parser.add_argument('--compare', type=Path)
    main(parser.parse_args())

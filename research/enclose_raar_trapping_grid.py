"""Attempt outward forward closure; report failures without attraction claims.

A closed run is a candidate certificate requiring an independent accounting
and source audit. The finite grid, exact window inputs, arithmetic assumptions,
positive image margin and all progress records must accompany such an audit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import mpmath as mp
import numpy as np
from check_raar_coverings_ad import exact
from propose_raar_trapping_bitset import mark, occupied_keys
from raar_interval_grid import advance


def indices(lower, upper, origin, h):
    a = np.floor(np.nextafter(np.nextafter(lower-origin, -np.inf)/h, -np.inf)).astype(np.int64)
    b = np.floor(np.nextafter(np.nextafter(upper-origin, np.inf)/h, np.inf)).astype(np.int64)
    return a, b


def mark_range(mask, a, b, shape):
    assert (a >= 0).all() and (b < shape).all()
    for i in range(a[0], b[0]+1):
        for j in range(a[1], b[1]+1):
            for k in range(a[2], b[2]+1):
                mark(mask, (i*shape[1]+j)*shape[2]+k)


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    begun = time.time()
    names = ['enclose_raar_trapping_grid.py', 'raar_interval_grid.py', 'raar_interval64.py',
             'raar_geometric_enclosure.py', 'propose_raar_trapping_bitset.py',
             'search_raar_circles.py', 'check_raar_coverings_ad.py']
    hashes = {}
    for name in names:
        path = Path(__file__).with_name(name)
        (args.out/name).write_bytes(path.read_bytes())
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    mp.iv.dps = 100
    reference = root/'results/raar-circles/render/case-1051-animation.npz'
    window_file = root/'results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json'
    d = np.load(reference)
    z, p, h = d['states'], d['params'], args.step
    assert 0 < p[0] < 1 and p[1] > 0 and h > 0
    assert 0 <= args.max_depth <= 10 and args.margin > 0
    assert args.reference_count >= 0 and args.checkpoint_every > 0 and args.generations > 0
    origin = np.floor((z.min(0)-args.domain_pad)/h)*h
    shape = np.ceil((z.max(0)+args.domain_pad-origin)/h).astype(np.int64)
    cells = int(np.prod(shape))
    assert cells <= args.max_grid_cells
    mask = np.zeros((cells+7)//8, np.uint8)
    window_boxes = []
    windows = json.loads(window_file.read_text())
    for c, B, r in zip(windows['centers'], windows['bases'], windows['radii']):
        bounds = []
        for i in range(3):
            radius = sum(abs(exact(B[i][j])*exact(r[j])) for j in range(3))
            lo, hi = exact(c[i])-radius, exact(c[i])+radius
            bounds.append([np.nextafter(float(lo.a), -np.inf), np.nextafter(float(hi.b), np.inf)])
            assert exact(bounds[-1][0]) < lo.a and exact(bounds[-1][1]) > hi.b
        box = np.array(bounds)
        a, b = indices(box[:, 0], box[:, 1], origin, h)
        mark_range(mask, a, b, shape)
        window_boxes.append({'box': bounds, 'lower_index': a.tolist(), 'upper_index': b.tolist()})
    reference_ids = np.linspace(0, len(z)-1, min(args.reference_count, len(z)), dtype=int) if args.reference_count else np.array([], dtype=int)
    for ident in reference_ids:
        a, b = indices(z[ident], z[ident], origin, h)
        mark_range(mask, a, b, shape)
    frontier = occupied_keys(mask)
    initial = len(frontier)
    trig = np.array([[np.nextafter(float(v.a), -np.inf), np.nextafter(float(v.b), np.inf)]
                     for v in [mp.iv.sin(exact(p[2])), mp.iv.cos(exact(p[2]))]])
    for interval, v in zip(trig, [mp.iv.sin(exact(p[2])), mp.iv.cos(exact(p[2]))]):
        assert exact(interval[0]) < v.a and exact(interval[1]) > v.b
    protocol = {'parameters': p.tolist(), 'trig_intervals': trig.tolist(), 'step': h,
                'origin': origin.tolist(), 'shape': shape.tolist(), 'grid_cells': cells,
                'initial_occupied': initial, 'reference_count': len(reference_ids),
                'curvature_tolerance': args.curvature_tolerance, 'max_depth': args.max_depth,
                'positive_image_margin': args.margin, 'window_boxes': window_boxes,
                'source_sha256': hashes, 'reference_sha256': hashlib.sha256(reference.read_bytes()).hexdigest(),
                'windows_sha256': hashlib.sha256(window_file.read_bytes()).hexdigest(),
                'arithmetic': 'Outward binary64 nextafter operations and correctly rounded platform sqrt; no fastmath. Exact binary64 parameter lifts; 100-digit interval trigonometric seeds.',
                'scope': 'An attempted finite forward closure. No completed invariant-set, attraction or novelty claim before a successful independent audit.'}
    (args.out/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n', encoding='utf-8')
    np.savez_compressed(args.out/'initial.npz', source_keys=frontier, reference_ids=reference_ids,
                        params=p, trig=trig, origin=origin, shape=shape, step=h, image_margin=args.margin)
    print(json.dumps({'stage': 'initial', 'cells': initial, 'seconds': time.time()-begun}), flush=True)
    history, flag, position, key = [], 5, 0, -1
    occupied_count = initial
    for generation in range(args.generations):
        new, flag, position, key, axis, splits, leaves, largest = advance(
            frontier, mask, shape, origin, h, p, trig, args.margin, args.max_new,
            args.curvature_tolerance, args.max_depth)
        occupied_count += len(new)
        row = {'generation': generation, 'sources': len(frontier), 'processed': int(position),
               'added': len(new), 'status': int(flag), 'axis_lower_bound': float(axis),
               'splits': int(splits), 'leaves': int(leaves), 'largest_remainder_over_step': float(largest),
               'occupied': occupied_count, 'seconds': time.time()-begun}
        history.append(row)
        (args.out/'progress.json').write_text(json.dumps(history, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(row), flush=True)
        if not flag:
            frontier = new
        if flag or not len(frontier) or (generation+1) % args.checkpoint_every == 0 or generation+1 == args.generations:
            np.savez_compressed(args.out/f'generation-{generation:04}.npz', packed_mask=mask,
                                pending=frontier, newly_added=new, next_generation=generation+1,
                                complete_generation=not bool(flag), processed_position=position)
        if flag or not len(frontier):
            break
    # A source snapshot changing during execution invalidates provenance.
    for name, digest in hashes.items():
        assert hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest, name
    report = {'status': 'closed_pending_audit' if flag == 0 and not len(frontier) else 'not_closed',
              'reason': {0: 'Closure or generation limit', 1: 'Regularity enclosure unresolved at depth limit',
                         2: 'Undefined or nonfinite center evaluation', 3: 'Proposed image bounds leave the declared domain; actual escape is not established',
                         4: 'New-cell budget exceeded', 5: 'No generations attempted'}[int(flag)],
              'initial_occupied': initial, 'final_occupied': occupied_count,
              'processed_cells': sum(v['processed'] for v in history), 'failing_key': int(key),
              'history': history, 'protocol_sha256': hashlib.sha256((args.out/'protocol.json').read_bytes()).hexdigest(),
              'scope': protocol['scope']}
    (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'history'}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--step', type=float, required=True)
    parser.add_argument('--reference-count', type=int, default=0)
    parser.add_argument('--domain-pad', type=float, default=.35)
    parser.add_argument('--max-grid-cells', type=int, default=4000000000)
    parser.add_argument('--generations', type=int, default=120)
    parser.add_argument('--max-new', type=int, default=12000000)
    parser.add_argument('--curvature-tolerance', type=float, default=.15)
    parser.add_argument('--max-depth', type=int, default=4)
    parser.add_argument('--margin', type=float, default=2.**-35)
    parser.add_argument('--checkpoint-every', type=int, default=5)
    main(parser.parse_args())

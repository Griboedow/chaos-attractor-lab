"""Resume a numerical proposal after enlarging positive coordinate bounds.

Keeps the exact origin and step. Replays the interrupted generation with an
empty mask to recover a superset of every newly added, unprocessed cell. This
is a numerical proposal continuation, never an interval certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from numba import njit
from propose_raar_trapping_bitset import mark, occupied_keys
from raar_adaptive_grid import advance


@njit(cache=True)
def rekey(keys, before, after):
    out = np.empty_like(keys)
    for n, key in enumerate(keys):
        i, j, k = key//(before[1]*before[2]), (key//before[2])%before[1], key%before[2]
        out[n] = (i*after[1]+j)*after[2]+k
    return out


@njit(cache=True)
def set_all(mask, keys):
    for key in keys:
        mark(mask, key)


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    begun = time.time()
    old = json.loads((args.input/'report.json').read_text())
    d = np.load(args.input/'proposal.npz')
    p, origin, before, h = [d[k] for k in ['params', 'origin', 'shape', 'step']]
    h = float(h)
    after = before+np.array(args.extend_upper, dtype=np.int64)
    assert (after >= before).all() and np.any(after > before)
    helper = Path(__file__).with_name('raar_adaptive_grid.py')
    assert hashlib.sha256(helper.read_bytes()).hexdigest() == old['adaptive_source_sha256']
    (args.out/helper.name).write_bytes(helper.read_bytes())
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    pending = d['pending']
    empty = np.zeros((int(np.prod(before))+7)//8, np.uint8)
    replay = advance(pending, empty, before, origin, h, p, old['image_buffer'], 12000000,
                     old['curvature_tolerance'], old['max_depth'])
    previous_last = old['history'][-1]
    assert int(replay[1]) == previous_last['status']
    assert int(replay[2]) == old['failing_position'] and int(replay[3]) == old['failing_key']
    frontier_old = np.unique(np.r_[pending[int(replay[2]):], replay[0]])
    occupied = d['occupied']
    # All replay targets were already inserted into the failed run's mask.
    assert np.isin(frontier_old, occupied, assume_unique=True).all()
    mask = np.zeros((int(np.prod(after))+7)//8, np.uint8)
    set_all(mask, rekey(occupied, before, after))
    frontier = rekey(frontier_old, before, after)
    protocol = {'input_report_sha256': hashlib.sha256((args.input/'report.json').read_bytes()).hexdigest(),
                'input_proposal_sha256': hashlib.sha256((args.input/'proposal.npz').read_bytes()).hexdigest(),
                'old_shape': before.tolist(), 'new_shape': after.tolist(), 'origin': origin.tolist(),
                'step': h, 'input_occupied': len(occupied), 'recovered_frontier': len(frontier),
                'replayed_new_keys': len(replay[0]), 'replay_failure_position': int(replay[2]),
                'helper_sha256': old['adaptive_source_sha256'],
                'scope': 'Continuation of a sampled image proposal, not an invariant-set proof.'}
    (args.out/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(protocol), flush=True)
    del occupied, empty, frontier_old, pending
    history, flag, position, key = [], 5, 0, -1
    for generation in range(args.generations):
        new, flag, position, key, axis, splits, leaves, largest = advance(
            frontier, mask, after, origin, h, p, old['image_buffer'], 12000000,
            old['curvature_tolerance'], old['max_depth'])
        row = {'generation': generation, 'sources': len(frontier), 'processed': int(position),
               'added': len(new), 'status': int(flag), 'minimum_sampled_axis': float(axis),
               'subdivision_splits': int(splits), 'subdivision_leaves': int(leaves),
               'largest_leaf_remainder_over_grid_step': float(largest), 'seconds': time.time()-begun}
        history.append(row)
        print(json.dumps(row), flush=True)
        (args.out/'progress.json').write_text(json.dumps(history, indent=2)+'\n', encoding='utf-8')
        if flag:
            break
        frontier = new
        if not len(frontier):
            break
    occupied = occupied_keys(mask)
    np.savez_compressed(args.out/'proposal.npz', params=p, origin=origin, shape=after, step=h,
                        occupied=occupied, pending=frontier, initial_reference_keys=rekey(d['initial_reference_keys'], before, after))
    report = dict(old)
    report.update(status='numerically_closed' if flag == 0 and not len(frontier) else 'not_closed',
                  continuation_protocol=protocol, history=history, final_occupied=len(occupied),
                  failing_key=int(key), failing_position=int(position),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    report['reason'] = {0: 'Sampled closure or generation limit', 1: 'Source cell intersects B axis',
                        2: 'Undefined sampled evaluation', 3: 'Proposed image bounds exceed the enlarged domain',
                        4: 'New-cell budget', 5: 'No generations'}[int(flag)]
    (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': report['status'], 'final_occupied': len(occupied), 'flag': int(flag)}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--extend-upper', type=int, nargs=3, default=[0, 0, 128])
    parser.add_argument('--generations', type=int, default=240)
    main(parser.parse_args())

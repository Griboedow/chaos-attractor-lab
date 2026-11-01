"""Validate one declared source-key range against a monotonically growing mask."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import numba
from raar_interval_offset_grid import advance


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    begun = time.time()
    protocol = json.loads((args.input/'protocol.json').read_text())
    files = ['audit_raar_trapping_batch.py', 'raar_interval_offset_grid.py', 'raar_geometric_scalar.py',
             'raar_interval_scalar.py', 'propose_raar_trapping_bitset.py', 'search_raar_circles.py']
    hashes = {}
    for name in files:
        path = Path(__file__).with_name(name)
        (args.out/name).write_bytes(path.read_bytes())
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    keys = np.load(args.input/'keys.npy', mmap_mode='r')
    mask = np.load(args.input/'mask.npy').copy()
    shape = np.array(protocol['shape'], dtype=np.int64)
    origin = np.array(protocol['origin'])
    p, trig = np.array(protocol['parameters']), np.array(protocol['trig_intervals'])
    offset = np.array(protocol['grid_offset'], dtype=np.int64)
    start, stop = args.start, args.stop if args.stop is not None else len(keys)
    assert 0 <= start < stop <= len(keys)
    selected = keys[start:stop]
    if args.samples:
        selected = np.asarray(selected[np.linspace(0, len(selected)-1, min(args.samples, len(selected)), dtype=int)])
    declaration = {'input_protocol_sha256': hashlib.sha256((args.input/'protocol.json').read_bytes()).hexdigest(),
                   'start': start, 'stop': stop, 'sample_count': args.samples, 'source_count': len(selected),
                   'max_depth': args.max_depth, 'source_sha256': hashes,
                   'runtime': {'python': sys.version, 'numpy': np.__version__, 'numba': numba.__version__},
                   'cache': 'Disabled for all outward arithmetic and raster kernels.',
                   'scope': 'One source batch only. A global closure claim requires a separate coverage, union and pending-source audit.'}
    (args.out/'declaration.json').write_text(json.dumps(declaration, indent=2)+'\n', encoding='utf-8')
    rows, parts = [], []
    print(json.dumps({'stage': 'start', 'selected': len(selected), 'start': start, 'stop': stop}), flush=True)
    flag = 0
    for local_start in range(0, len(selected), args.chunk_size):
        source = np.asarray(selected[local_start:local_start+args.chunk_size])
        result = advance(source, mask, shape, origin, protocol['step'], p, trig,
                         protocol['positive_image_margin'], args.max_new, args.max_depth, offset)
        new, flag, processed, key, axis, splits, leaves, depth, bbox = result
        path = args.out/f'new-{local_start:09}.npy'
        np.save(path, new)
        row = {'local_start': local_start, 'sources': len(source), 'processed': int(processed),
               'source_keys_sha256': hashlib.sha256(source.tobytes()).hexdigest(),
               'new_file': path.name, 'new_file_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'added': len(new), 'status': int(flag), 'failing_key': int(key),
               'axis_lower_bound': float(axis), 'splits': int(splits), 'leaves': int(leaves),
               'maximum_accepted_depth': int(depth), 'bounding_box_acceptances': int(bbox),
               'seconds': time.time()-begun}
        rows.append(row)
        (args.out/'progress.json').write_text(json.dumps(rows, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(row), flush=True)
        if flag:
            break
    for name, digest in hashes.items():
        assert hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest, name
    report = {'completed': flag == 0 and sum(v['processed'] for v in rows) == len(selected),
              'processed': sum(v['processed'] for v in rows), 'added': sum(v['added'] for v in rows),
              'axis_lower_bound': min(v['axis_lower_bound'] for v in rows), 'history': rows,
              'declaration_sha256': hashlib.sha256((args.out/'declaration.json').read_bytes()).hexdigest(),
              'scope': declaration['scope']}
    (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'history'}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--stop', type=int)
    parser.add_argument('--samples', type=int, default=0)
    parser.add_argument('--chunk-size', type=int, default=100000)
    parser.add_argument('--max-depth', type=int, default=6)
    parser.add_argument('--max-new', type=int, default=1000000)
    main(parser.parse_args())

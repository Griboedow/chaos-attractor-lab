"""Independent input, mask, window-inclusion and exhaustive-range accounting.

This auditor does not replace inspection or rerunning of the interval kernels.
It independently checks their immutable inputs, exact rational window coverage,
all source-range hashes, subdivision accounting and every returned new-cell file.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import mpmath as mp
import numpy as np


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(part)
    return digest.hexdigest()


def exact(x):
    value = Q(float(x))
    return mp.iv.mpf(value.numerator)/value.denominator


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    protocol_file = args.input/'protocol.json'
    protocol = json.loads(protocol_file.read_text())
    assert sha(args.input/'mask.npy') == protocol['mask_sha256']
    assert sha(args.input/'keys.npy') == protocol['keys_sha256']
    mask = np.load(args.input/'mask.npy', mmap_mode='r')
    keys = np.load(args.input/'keys.npy', mmap_mode='r')
    assert mask.dtype == np.uint8 and keys.dtype == np.int64 and mask.ndim == keys.ndim == 1
    shape = np.array(protocol['shape'], dtype=np.int64)
    count = int(np.prod(shape))
    assert len(mask) == (count+7)//8
    pointer = 0
    # Independent NumPy unpacking, rather than the compiled packed-grid helpers.
    for begin in range(0, len(mask), 1_000_000):
        bits = np.unpackbits(np.asarray(mask[begin:begin+1_000_000]), bitorder='little')
        locations = np.flatnonzero(bits)+8*begin
        assert np.array_equal(locations, keys[pointer:pointer+len(locations)])
        assert not len(locations) or locations[-1] < count
        pointer += len(locations)
    assert pointer == len(keys) == protocol['initial_cells']
    window_file = root/'results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json'
    assert sha(window_file) == protocol['windows_sha256']
    windows = json.loads(window_file.read_text())
    assert protocol['parameters'] == windows['params']
    assert len(windows['centers']) == len(windows['bases']) == len(windows['radii']) == 7
    origin, h = [Q(float(v)) for v in protocol['origin']], Q(float(protocol['step']))
    offset = protocol['grid_offset']
    assert h > 0 and 0 < protocol['parameters'][0] < 1 and protocol['parameters'][1] > 0
    margin = Q(float(protocol['positive_image_margin']))
    assert margin > 0
    exact_window_cells = []
    for center, basis, radii in zip(windows['centers'], windows['bases'], windows['radii']):
        a, b = [], []
        for i in range(3):
            radius = sum(abs(Q(float(basis[i][j]))*Q(float(radii[j]))) for j in range(3))
            lower = (Q(float(center[i]))-radius-origin[i])/h
            upper = (Q(float(center[i]))+radius-origin[i])/h
            a.append(lower.numerator//lower.denominator+offset[i])
            b.append(upper.numerator//upper.denominator+offset[i])
        checked = 0
        for i in range(a[0], b[0]+1):
            for j in range(a[1], b[1]+1):
                for k in range(a[2], b[2]+1):
                    assert 0 <= i < shape[0] and 0 <= j < shape[1] and 0 <= k < shape[2]
                    key = int((i*shape[1]+j)*shape[2]+k)
                    assert (int(mask[key//8]) >> (key%8)) & 1
                    checked += 1
        exact_window_cells.append({'lower_index': a, 'upper_index': b, 'checked_cells': checked})
    mp.iv.dps = 120
    theta = exact(protocol['parameters'][2])
    for saved, actual in zip(protocol['trig_intervals'], [mp.iv.sin(theta), mp.iv.cos(theta)]):
        assert exact(saved[0]) < actual.a and exact(saved[1]) > actual.b
    # The earlier located periodic roots now lie in the trapping candidate.
    # Recheck their tiny Cartesian boxes directly with exact rational indexing.
    cycle_membership = []
    cycle_paths = ['results/raar-cycle-structure/cycle-proof/base-unstable-period-7.json',
                   'results/raar-symbolic-dynamics/runs/raar-symbolic-cycle-proof01/symbolic-period-49.json']
    for relative in cycle_paths:
        path = root/relative
        certificate = json.loads(path.read_text())
        assert certificate['parameters'] == protocol['parameters']
        cube_bounds = []
        for i, entry in enumerate(certificate['box']):
            ends = [Q(v.strip()) for v in entry.strip('[]').split(',')]
            raw = [(v-origin[i])/h for v in ends]
            cube_bounds.append([v.numerator//v.denominator+offset[i] for v in raw])
        cells = 0
        for i in range(cube_bounds[0][0], cube_bounds[0][1]+1):
            for j in range(cube_bounds[1][0], cube_bounds[1][1]+1):
                for k in range(cube_bounds[2][0], cube_bounds[2][1]+1):
                    assert 0 <= i < shape[0] and 0 <= j < shape[1] and 0 <= k < shape[2]
                    key = int((i*shape[1]+j)*shape[2]+k)
                    assert (int(mask[key//8]) >> (key%8)) & 1
                    cells += 1
        row = {'period': certificate['period'], 'certificate': relative,
               'certificate_sha256': sha(path), 'cube_bounds': cube_bounds, 'checked_cells': cells}
        if certificate['period'] == 7:
            matrix = [[mp.iv.mpf([v.strip() for v in entry.strip('[]').split(',')]) for entry in line]
                      for line in certificate['monodromy']]
            a, b, c = matrix
            determinant = a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
            assert determinant.b < 0
            row['full_period_determinant_interval'] = str(determinant)
        cycle_membership.append(row)
    batches, all_new, shared_source_hashes = [], [], None
    total_splits = total_leaves = total_bbox = 0
    minimum_axis = float('inf')
    for directory in args.batches:
        declaration_path, report_path = directory/'declaration.json', directory/'report.json'
        declaration, report = json.loads(declaration_path.read_text()), json.loads(report_path.read_text())
        assert report['completed'] and report['declaration_sha256'] == sha(declaration_path)
        assert declaration['sample_count'] == 0
        assert declaration['input_protocol_sha256'] == sha(protocol_file)
        assert declaration['input_mask_sha256'] == protocol['mask_sha256']
        assert declaration['input_keys_sha256'] == protocol['keys_sha256']
        if shared_source_hashes is None:
            shared_source_hashes = declaration['source_sha256']
        assert declaration['source_sha256'] == shared_source_hashes
        for name, digest in shared_source_hashes.items():
            assert sha(directory/name) == digest, (directory, name)
        start, stop = declaration['start'], declaration['stop']
        assert 0 <= start < stop <= len(keys)
        assert declaration['source_count'] == stop-start
        at, added, batch_new = 0, 0, []
        for row in report['history']:
            assert row['status'] == 0 and row['processed'] == row['sources'] > 0
            assert row['local_start'] == at and start+at+row['sources'] <= stop
            source = np.asarray(keys[start+at:start+at+row['sources']])
            assert hashlib.sha256(source.tobytes()).hexdigest() == row['source_keys_sha256']
            assert row['leaves'] == row['sources']+7*row['splits']
            assert 0 <= row['bounding_box_acceptances'] <= row['leaves']
            assert 0 <= row['maximum_accepted_depth'] <= declaration['max_depth']
            assert np.isfinite(row['axis_lower_bound']) and row['axis_lower_bound'] > 0
            minimum_axis = min(minimum_axis, row['axis_lower_bound'])
            total_splits += row['splits']; total_leaves += row['leaves']; total_bbox += row['bounding_box_acceptances']
            path = directory/row['new_file']
            assert path.resolve().parent == directory.resolve()
            assert sha(path) == row['new_file_sha256']
            new = np.load(path)
            assert new.dtype == np.int64 and new.ndim == 1 and len(new) == row['added']
            if len(new):
                assert new.min() >= 0 and new.max() < count
                assert np.all(((np.asarray(mask[new//8]) >> (new%8)) & 1) == 0)
                batch_new.append(new)
            at += row['sources']; added += len(new)
        assert at == stop-start == report['processed'] and added == report['added']
        if batch_new:
            merged = np.concatenate(batch_new)
            assert len(np.unique(merged)) == len(merged)
            all_new.append(merged)
        batches.append({'start': start, 'stop': stop, 'processed': at, 'added': added,
                        'declaration_sha256': sha(declaration_path), 'report_sha256': sha(report_path)})
    batches.sort(key=lambda v: v['start'])
    at = 0
    for batch in batches:
        assert batch['start'] == at
        at = batch['stop']
    assert at == len(keys)
    new_keys = np.unique(np.concatenate(all_new)) if all_new else np.array([], dtype=np.int64)
    np.save(args.out/'pending-new-keys.npy', new_keys)
    volume = len(keys)*h**3
    result = {'status': 'closed_accounting_verified' if not len(new_keys) else 'additional_cells_require_processing',
              'source_cells': len(keys), 'new_unprocessed_cells': len(new_keys), 'batches': batches,
              'accepted_leaves': total_leaves, 'subdivision_splits': total_splits,
              'bounding_box_acceptances': total_bbox, 'axis_lower_bound': minimum_axis,
              'image_margin_exact': [margin.numerator, margin.denominator],
              'candidate_volume_exact': [volume.numerator, volume.denominator],
              'candidate_volume_approximate': float(volume), 'exact_window_inclusion': exact_window_cells,
              'located_cycle_inclusion': cycle_membership,
              'protocol_sha256': sha(protocol_file), 'interval_source_sha256': shared_source_hashes,
              'auditor_sha256': sha(Path(__file__)),
              'scope': 'Input and exhaustive computational coverage audit. Invariance relies on the separately inspected outward geometric enclosure and raster kernels; attraction and entropy are mathematical consequences under their stated arithmetic assumptions. No transitivity or historical-priority claim.'}
    (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps({k: v for k, v in result.items() if k not in ['batches', 'interval_source_sha256', 'exact_window_inclusion']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--batches', type=Path, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    main(parser.parse_args())

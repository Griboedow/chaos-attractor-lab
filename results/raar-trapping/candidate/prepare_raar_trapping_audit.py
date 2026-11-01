"""Freeze an explicit cubical candidate and certify inclusion of graph windows."""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from check_raar_coverings_ad import exact
from enclose_raar_trapping_grid import indices
from propose_raar_trapping_bitset import mark, marked, occupied_keys
from resume_raar_adaptive_grid import set_all


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(part)
    return digest.hexdigest()


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    report = json.loads((args.input/'report.json').read_text())
    data = np.load(args.input/'proposal.npz')
    p, origin, shape, h = [data[k] for k in ['params', 'origin', 'shape', 'step']]
    h = float(h)
    offset = np.array(report.get('grid_offset', [0, 0, 0]), dtype=np.int64)
    original_keys = data['occupied']
    assert np.all(np.diff(original_keys) > 0)
    mask = np.zeros((int(np.prod(shape))+7)//8, np.uint8)
    set_all(mask, original_keys)
    mp.iv.dps = 100
    window_source = root/'results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json'
    windows = json.loads(window_source.read_text())
    boxes, additions = [], []
    for c, B, r in zip(windows['centers'], windows['bases'], windows['radii']):
        box = np.empty((3, 2))
        for i in range(3):
            radius = sum(abs(exact(B[i][j])*exact(r[j])) for j in range(3))
            lo, hi = exact(c[i])-radius, exact(c[i])+radius
            box[i] = [np.nextafter(float(lo.a), -np.inf), np.nextafter(float(hi.b), np.inf)]
            assert exact(box[i, 0]) < lo.a and exact(box[i, 1]) > hi.b
        a, b = indices(box[:, 0], box[:, 1], origin, h)
        a, b = a+offset, b+offset
        assert (a >= 0).all() and (b < shape).all()
        for i in range(a[0], b[0]+1):
            for j in range(a[1], b[1]+1):
                for k in range(a[2], b[2]+1):
                    key = (i*shape[1]+j)*shape[2]+k
                    if not marked(mask, key):
                        mark(mask, key); additions.append(int(key))
        boxes.append({'bounds': box.tolist(), 'lower_storage_index': a.tolist(), 'upper_storage_index': b.tolist()})
    keys = occupied_keys(mask)
    assert len(keys) == len(original_keys)+len(additions)
    trig = np.array([[np.nextafter(float(v.a), -np.inf), np.nextafter(float(v.b), np.inf)]
                     for v in [mp.iv.sin(exact(p[2])), mp.iv.cos(exact(p[2]))]])
    for interval, v in zip(trig, [mp.iv.sin(exact(p[2])), mp.iv.cos(exact(p[2]))]):
        assert exact(interval[0]) < v.a and exact(interval[1]) > v.b
    np.save(args.out/'mask.npy', mask)
    np.save(args.out/'keys.npy', keys)
    protocol = {'parameters': p.tolist(), 'origin': origin.tolist(), 'shape': shape.tolist(),
                'step': h, 'grid_offset': offset.tolist(), 'trig_intervals': trig.tolist(),
                'positive_image_margin': 2.**-35, 'initial_cells': len(keys),
                'window_cover': boxes, 'added_for_exact_window_inclusion': additions,
                'input_report_sha256': sha(args.input/'report.json'),
                'input_proposal_sha256': sha(args.input/'proposal.npz'),
                'windows_sha256': sha(window_source), 'mask_sha256': sha(args.out/'mask.npy'),
                'keys_sha256': sha(args.out/'keys.npy'),
                'scope': 'An explicit candidate set including the seven certified graph windows. No forward-invariance or attraction assertion yet.'}
    (args.out/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n', encoding='utf-8')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps({k: v for k, v in protocol.items() if k not in ['window_cover', 'parameters', 'trig_intervals']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    main(parser.parse_args())

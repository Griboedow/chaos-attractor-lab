"""Outward arbitrary-precision covering relations for G=F^4.

Every closed source cube and both exit faces are covered by a complete binary
subdivision tree. This certifies a bounded symbolic subsystem, not attraction.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from product_binary_interval import exact, parameters, inverse, mm, mv, fourth


def main(args):
    sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    mp.mp.dps = mp.iv.dps = args.digits
    iv = mp.iv
    begin = time.time()
    args.out.mkdir(parents=True, exist_ok=False)
    save = lambda p, v: p.write_text(json.dumps(v, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    sources = {}
    for name in [Path(__file__).name, 'product_binary_interval.py']:
        dest = args.out/name
        dest.write_bytes(Path(__file__).with_name(name).read_bytes())
        sources[name] = sha(dest)
    with np.load(args.input) as data:
        raw = {k: data[k] for k in ['params', 'points', 'centers', 'bases', 'radii', 'edges', 'orientation']}
    par = parameters(raw['params'], raw['points'], iv)
    centers = [[exact(v, iv) for v in row] for row in raw['centers']]
    matrices, inverses, exact_windows = [], [], []
    for b, r in zip(raw['bases'], raw['radii']):
        mat = [[Fraction(float(b[i, j]))*Fraction(float(r[j])) for j in range(3)] for i in range(3)]
        inv = inverse(mat)
        assert mm(mat, inv) == [[Fraction(i == j) for j in range(3)] for i in range(3)]
        cv = lambda a: [[iv.mpf(v.numerator)/iv.mpf(v.denominator) for v in row] for row in a]
        matrices.append(cv(mat))
        inverses.append(cv(inv))
        exact_windows.append(dict(matrix=[[str(v) for v in row] for row in mat], inverse=[[str(v) for v in row] for row in inv]))
    unit = [iv.mpf([-1, 1])]*3
    sep = mv(inverses[0], [centers[1][j]-centers[0][j] for j in range(3)])[0]+mv(mm(inverses[0], matrices[1]), unit)[0]
    assert sep.b < -1 or sep.a > 1
    save(args.out/'input.json', dict(**{k: v.tolist() for k, v in raw.items()}, exact_windows=exact_windows, proposal_sha256=sha(args.input), window01_separation=str(sep)))

    def enclose(i, j, mid, rad):
        q = [iv.mpf(v) for v in mid]
        delta = [iv.mpf([-v, v]) for v in rad]
        off, spread = mv(matrices[i], q), mv(matrices[i], delta)
        center = [centers[i][k]+off[k] for k in range(3)]
        box = [center[k]+spread[k] for k in range(3)]
        _, jac, radii = fourth(box, par, iv)
        assert all(r.a > 0 for r in radii)
        image, _, _ = fourth(center, par, iv)
        out = mv(inverses[j], [image[k]-centers[j][k] for k in range(3)])
        error = mv(mm(mm(inverses[j], jac), matrices[i]), delta)
        return [out[k]+error[k] for k in range(3)], min(mp.mpf(r.a) for r in radii)

    def accepted(y, kind, orientation):
        if kind == 'volume':
            return y[0].b < -1 or y[0].a > 1 or all(v.a > -1 and v.b < 1 for v in y[1:])
        sign = -orientation if kind == 'left' else orientation
        return (sign*y[0]).a > 1

    records = []
    for edge_id, (i0, j0) in enumerate(raw['edges']):
        i, j = int(i0), int(j0)
        orientation = int(raw['orientation'][edge_id])
        assert orientation in [-1, 1]
        scales = np.sum(abs(raw['bases'][i]*raw['radii'][i]), axis=0)
        for kind in ['volume', 'left', 'right']:
            mid, rad = [mp.mpf(0)]*3, [mp.mpf(1)]*3
            if kind != 'volume':
                mid[0] = mp.mpf(-1 if kind == 'left' else 1)
                rad[0] = mp.mpf(0)
            stack = [(mid, rad, '')]
            examined = leaves = depth = 0
            minimum = mp.inf
            last = time.time()
            path = args.out/f'edge-{edge_id}-{kind}.jsonl.gz'
            with gzip.open(path, 'wt', encoding='utf-8', newline='\n') as handle:
                while stack and examined < args.max_cells:
                    mid, rad, address = stack.pop()
                    examined += 1
                    depth = max(depth, len(address)//2)
                    y, axis = enclose(i, j, mid, rad)
                    if accepted(y, kind, orientation):
                        leaves += 1
                        minimum = min(minimum, axis)
                        handle.write(json.dumps(dict(path=address, image=[str(v) for v in y], minimum_radius=mp.nstr(axis, args.digits)), separators=(',', ':'))+'\n')
                    else:
                        split = max(range(3), key=lambda k: scales[k]*float(rad[k]))
                        assert rad[split] > mp.mpf(2)**-45, (edge_id, kind, address)
                        child = rad.copy()
                        child[split] /= 2
                        for sign in [-1, 1]:
                            cm = mid.copy()
                            cm[split] += sign*child[split]
                            stack.append((cm, child.copy(), address+str(split)+('0' if sign < 0 else '1')))
                    if time.time()-last > 15:
                        print(json.dumps(dict(edge=edge_id, kind=kind, examined=examined, pending=len(stack), seconds=time.time()-begin)), flush=True)
                        last = time.time()
            row = dict(edge=edge_id, source=i, target=j, kind=kind, status='passed' if not stack else 'incomplete', examined=examined, leaves=leaves, pending=len(stack), maximum_depth=depth, minimum_radius=mp.nstr(minimum, args.digits), certificate_sha256=sha(path))
            records.append(row)
            print(json.dumps(row), flush=True)
            save(args.out/'progress.json', records)
            if stack:
                save(args.out/f'edge-{edge_id}-{kind}-pending.json', [[list(map(str, m)), list(map(str, r)), p] for m, r, p in stack])
                break
        if stack:
            break
    passed = len(records) == 3*len(raw['edges']) and all(r['status'] == 'passed' for r in records)
    result = dict(status='passed' if passed else 'incomplete', precision_digits=args.digits, mpmath_version=mp.__version__, input_sha256=sha(args.out/'input.json'), sources=sources, checks=records, seconds=time.time()-begin, scope='Outward covering relations for G=F^4 with disjoint N0,N1. No trapping region or attraction certificate. Symbolic consequences require a complete coverage and source audit. No historical priority claim.')
    save(args.out/'report.json', result)
    print(json.dumps(dict(status=result['status'], seconds=result['seconds'])), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['input', 'out', 'mpmath-path']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--digits', type=int, default=80)
    p.add_argument('--max-cells', type=int, default=2000)
    main(p.parse_args())

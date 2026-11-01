"""Certify the sign of det D(F^4) on every covering window.

Multiply four one-step determinants instead of taking the determinant of an
interval matrix product. Both are valid, but the former has less dependency.
The full-site engine shares interval arithmetic and windows, not the formula
or derivative implementation, with the binary-factorized engine.
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
from product_binary_interval import parameters, model, mv
from certify_product_cycles import det


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    if args.mpmath_path:
        sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    import verify_product_coverings as full
    full.mp = mp
    mp.mp.dps = mp.iv.dps = args.digits
    iv = mp.iv
    exact = lambda v: iv.mpf(Fraction(v).numerator)/Fraction(v).denominator
    enc = lambda v: dict(decimal=str(v), binary=[list(t) for t in v._mpi_])
    args.out.mkdir(parents=True, exist_ok=False)
    raw = json.loads((args.source/'input.json').read_text(encoding='utf-8'))
    original = json.loads((args.source/'report.json').read_text(encoding='utf-8'))
    assert original['status'] == 'passed' and original['input_sha256'] == sha(args.source/'input.json')
    assert raw['edges'] == [[0, 0], [0, 1]]+[[i, (i+1) % 23] for i in range(1, 23)]
    par = parameters(raw['params'], np.array(raw['points']), iv)
    independent = full.Model(raw)
    step = independent.value_jacobian if args.engine == 'full' else lambda x: model(x, par, iv)
    start, rows, total = time.time(), [], iv.mpf(1)
    for i, window in enumerate(raw['exact_windows']):
        center = list(map(exact, raw['centers'][i]))
        S = [[exact(v) for v in row] for row in window['matrix']]
        stack, leaves, lower, upper, minimum = [''], 0, None, None, None
        path = args.out/f'window-{i}.jsonl.gz'
        with gzip.open(path, 'wt', encoding='utf-8', newline='\n') as handle:
            while stack:
                address = stack.pop()
                mid, rad = full.cell(address, 'volume')
                local = [iv.mpf([(exact(a)-exact(b)).a, (exact(a)+exact(b)).b]) for a, b in zip(mid, rad)]
                spread = mv(S, local)
                box = [a+b for a, b in zip(center, spread)]
                product, radii, one_step = iv.mpf(1), [], []
                for _ in range(4):
                    natural, jac, radius = step(box)
                    assert radius.a > 0
                    value = det(jac)
                    product *= value
                    radii.append(radius)
                    one_step.append(value)
                    if args.engine == 'full':
                        middle = [iv.mpf((mp.mpf(v.a)+mp.mpf(v.b))/2) for v in box]
                        assert all(c.a >= v.a and c.b <= v.b for c, v in zip(middle, box))
                        base, _ = independent.field(middle)
                        change = mv(jac, [v-c for v, c in zip(box, middle)])
                        mean = [a+b for a, b in zip(base, change)]
                        box = [iv.mpf([max(a.a, b.a), min(a.b, b.b)]) for a, b in zip(natural, mean)]
                    else:
                        box = natural
                expected = -1 if i in [4, 5, 6] else 1
                passed = product.b < 0 if expected == -1 else product.a > 0
                if passed:
                    leaves += 1
                    lower = product.a if lower is None else min(lower, product.a)
                    upper = product.b if upper is None else max(upper, product.b)
                    minimum = min(v.a for v in radii) if minimum is None else min(minimum, *(v.a for v in radii))
                    handle.write(json.dumps(dict(path=address, determinant=enc(product), one_step=[enc(v) for v in one_step], radii=[enc(v) for v in radii]), separators=(',', ':'))+'\n')
                else:
                    if len(address)//2 >= 30 or leaves+len(stack) > args.max_cells:
                        raise ArithmeticError(f'Orientation subdivision budget exhausted at window {i}')
                    widths = [float(rad[j])*sum(abs(float(Fraction(v[j]))) for v in window['matrix']) for j in range(3)]
                    axis = max(range(3), key=lambda j: widths[j])
                    stack.extend([address+str(axis)+'0', address+str(axis)+'1'])
        interval = iv.mpf([lower, upper])
        total *= interval
        row = dict(window=i, leaves=leaves, determinant=enc(interval), minimum_radius=enc(minimum), sign=expected, file=path.name, sha256=sha(path))
        rows.append(row)
        print(json.dumps(dict(window=i, leaves=leaves, sign=expected, determinant=[float(lower), float(upper)])), flush=True)
    assert total.b < 0
    sources = {}
    for name in [Path(__file__).name, 'certify_product_cycles.py', 'product_binary_interval.py', 'verify_product_coverings.py']:
        p = Path(__file__).with_name(name)
        (args.out/name).write_bytes(p.read_bytes())
        sources[name] = sha(p)
    (args.out/'input.json').write_bytes((args.source/'input.json').read_bytes())
    report = dict(status='passed', engine=args.engine, digits=args.digits, mpmath=mp.__version__,
                  input_sha256=sha(args.source/'input.json'), covering_report_sha256=sha(args.source/'report.json'),
                  sources=sources, windows=rows, excursion_determinant=enc(total), seconds=time.time()-start,
                  scope='Every window is covered by a complete binary partition with a strict determinant sign. The 23-arrow covering excursion has a periodic realization with det D(F^92)<0. Least period 92 follows from N0/N1 separation and the even-power obstruction for periods 23 or 46. No attraction or intrinsic topological-conjugacy conclusion.')
    (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--mpmath-path', type=Path)
    p.add_argument('--digits', type=int, default=100)
    p.add_argument('--engine', choices=['factorized', 'full'], default='factorized')
    p.add_argument('--max-cells', type=int, default=20000)
    main(p.parse_args())

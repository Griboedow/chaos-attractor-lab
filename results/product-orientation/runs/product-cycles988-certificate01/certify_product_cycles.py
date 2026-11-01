"""Interval existence, least period and multipliers for binary product cycles.

The full-site mode recomputes with eight Gaussian weights, nested reflections
and forward automatic derivatives. It shares the proof driver, proposed center,
preconditioner and mpmath arithmetic with the factorized calculation.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from product_binary_interval import parameters, model, mm, mv


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def eye(ctx):
    return [[ctx.mpf(int(i == j)) for j in range(3)] for i in range(3)]


def det(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def iterate(x, period, ctx, step):
    total, orbit, radii = eye(ctx), [], []
    for _ in range(period):
        orbit.append(x)
        x, jac, r = step(x)
        assert r > 0
        total = mm(jac, total)
        radii.append(r)
    return x, total, orbit, radii


def main(args):
    if args.mpmath_path:
        sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    mp.mp.dps = mp.iv.dps = args.digits
    iv = mp.iv
    args.out.mkdir(parents=True, exist_ok=False)
    start, rows = time.time(), []
    raw = json.loads(args.windows.read_text(encoding='utf-8'))
    par = parameters(raw['params'], np.array(raw['points']), mp.mp)
    pari = parameters(raw['params'], np.array(raw['points']), iv)
    if args.engine == 'full':
        import verify_product_coverings as full
        full.mp = mp
        independent = full.Model(raw)
        step = independent.value_jacobian
        primary = json.loads((args.primary/'report.json').read_text(encoding='utf-8'))
        assert primary['status'] == 'passed'
        assert primary['windows_sha256'] == sha(args.windows)
    else:
        step = lambda x: model(x, pari, iv)
    decimal = lambda v: mp.nstr(v, args.digits+10)
    binary = lambda v: [list(t) for t in v._mpi_]
    enc = lambda v: dict(decimal=str(v), binary=binary(v))
    for ident in args.ids:
        source = args.source/f'cycle-{ident}.npz'
        with np.load(source) as data:
            assert np.array_equal(data['params'], raw['params'])
            assert np.array_equal(data['points'], raw['points'])
            seed, period = data['orbit'][0], len(data['orbit'])
        I = eye(mp.mp)
        if args.engine == 'factorized':
            x = [mp.mpf(float(v)) for v in seed]
            for iteration in range(30):
                y, M, _, _ = iterate(x, period, mp.mp, lambda z: model(z, par, mp.mp))
                residual = [a-b for a, b in zip(y, x)]
                if max(abs(v) for v in residual) < mp.mpf(10)**(-args.digits+30):
                    break
                delta = mp.lu_solve(mp.matrix([[M[i][j]-I[i][j] for j in range(3)] for i in range(3)]), mp.matrix(residual))
                x = [a-b for a, b in zip(x, delta)]
            else:
                raise ArithmeticError('Newton refinement failed')
            preconditioner = mp.inverse(mp.matrix([[M[i][j]-I[i][j] for j in range(3)] for i in range(3)]))
            center_text = [decimal(v) for v in x]
            inverse_text = [[decimal(preconditioner[i, j]) for j in range(3)] for i in range(3)]
            eig = sorted(mp.eig(mp.matrix(M), left=False, right=False), key=lambda v: mp.re(v))
            assert all(abs(mp.im(v)) < mp.mpf('1e-100') for v in eig)
            brackets = [[decimal(mp.re(v)-mp.mpf('1e-25')), decimal(mp.re(v)+mp.mpf('1e-25'))] for v in eig]
        else:
            prior = next(row for row in primary['cycles'] if row['id'] == ident)
            assert prior['source_sha256'] == sha(source) and prior['period'] == period
            center_text, inverse_text = prior['center_decimal'], prior['preconditioner_decimal']
            brackets = [v['endpoints_decimal'] for v in prior['multipliers']]
        center = [iv.mpf(v) for v in center_text]
        C = [[iv.mpf(v) for v in row] for row in inverse_text]
        radius = iv.mpf(args.radius)
        box = [iv.mpf([(v-radius).a, (v+radius).b]) for v in center]
        y0, _, _, _ = iterate(center, period, iv, step)
        _, M, orbit, radii = iterate(box, period, iv, step)
        II = eye(iv)
        H = [[M[i][j]-II[i][j] for j in range(3)] for i in range(3)]
        CH = mm(C, H)
        R = [[II[i][j]-CH[i][j] for j in range(3)] for i in range(3)]
        correction = mv(C, [a-b for a, b in zip(y0, center)])
        remainder = mv(R, [a-b for a, b in zip(box, center)])
        K = [a-b+c for a, b, c in zip(center, correction, remainder)]
        assert all(k.a > b.a and k.b < b.b for k, b in zip(K, box))
        contraction = max(sum(abs(v).b for v in row) for row in R)
        assert contraction < 1 and (det(C).a > 0 or det(C).b < 0)
        exclusions = []
        for d in range(1, period):
            if period % d:
                continue
            difference = [a-b for a, b in zip(orbit[d], box)]
            assert any(v.a > 0 or v.b < 0 for v in difference)
            exclusions.append(dict(divisor=d, differences=[enc(v) for v in difference]))
        coefficients = [-sum(M[i][i] for i in range(3)),
                        sum(M[i][i]*M[j][j]-M[i][j]*M[j][i] for i in range(3) for j in range(i+1, 3)), -det(M)]
        def polynomial(v):
            return ((v+coefficients[0])*v+coefficients[1])*v+coefficients[2]
        multipliers = []
        for lo, hi in brackets:
            left, right = polynomial(iv.mpf(lo)), polynomial(iv.mpf(hi))
            assert (left.b < 0 and right.a > 0) or (left.a > 0 and right.b < 0)
            enclosure = iv.mpf([iv.mpf(lo).a, iv.mpf(hi).b])
            stable = bool(enclosure.a > -1 and enclosure.b < 1)
            unstable = bool(enclosure.b < -1 or enclosure.a > 1)
            assert stable or unstable
            multipliers.append(dict(endpoints_decimal=[lo, hi], enclosure=enc(enclosure), at_left=enc(left), at_right=enc(right), stable=stable, unstable=unstable))
        assert all(iv.mpf(brackets[i][1]).b < iv.mpf(brackets[i+1][0]).a for i in range(2))
        local = mv([[full.exact(v) if args.engine == 'full' else iv.mpf(int(v.split('/')[0]))/iv.mpf(int(v.split('/')[1])) if '/' in v else iv.mpf(v) for v in row] for row in raw['exact_windows'][0]['inverse']], [box[j]-iv.mpf(float(raw['centers'][0][j])) for j in range(3)])
        in_window = all(v.a > -1 and v.b < 1 for v in local)
        row = dict(id=ident, period=period, source_sha256=sha(source), center_decimal=center_text,
                   preconditioner_decimal=inverse_text, radius=args.radius, box=[enc(v) for v in box],
                   krawczyk=[enc(v) for v in K], contraction=enc(contraction), preconditioner_determinant=enc(det(C)),
                   monodromy=[[enc(v) for v in r] for r in M], determinant=enc(det(M)),
                   characteristic_coefficients=[enc(v) for v in coefficients], characteristic_at_one=enc(polynomial(iv.mpf(1))),
                   multipliers=multipliers, proper_divisor_exclusions=exclusions,
                   orbit=[[enc(v) for v in r] for r in orbit], radii=[enc(v) for v in radii],
                   in_window_zero=bool(in_window), window_zero_coordinates=[enc(v) for v in local])
        rows.append(row)
        print(json.dumps(dict(id=ident, period=period, determinant=str(det(M)), in_window_zero=bool(in_window))), flush=True)
    sources = {}
    for name in [Path(__file__).name, 'product_binary_interval.py', 'verify_product_coverings.py']:
        path = Path(__file__).with_name(name)
        (args.out/name).write_bytes(path.read_bytes())
        sources[name] = sha(path)
    (args.out/'windows.json').write_bytes(args.windows.read_bytes())
    result = dict(status='passed', engine=args.engine, digits=args.digits, mpmath=mp.__version__,
                  windows_sha256=sha(args.windows), primary_sha256=sha(args.primary/'report.json') if args.primary else None,
                  sources=sources, cycles=rows, seconds=time.time()-start,
                  scope='Existence and uniqueness in each box, least period, three real hyperbolic multipliers and return determinant. Window-zero membership is tested separately. No attraction or equality with a numerical cloud; no general historical uniqueness conclusion.')
    (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['source', 'windows', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--primary', type=Path)
    p.add_argument('--mpmath-path', type=Path)
    p.add_argument('--ids', type=int, nargs='+', default=[14, 16, 23])
    p.add_argument('--engine', choices=['factorized', 'full'], default='factorized')
    p.add_argument('--digits', type=int, default=180)
    p.add_argument('--radius', default='1e-60')
    main(p.parse_args())

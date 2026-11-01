"""Independent arithmetic and map checks for the binary64 interval extension."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

import mpmath as mp
import numpy as np
import raar_interval64 as fast
import raar_geometric_enclosure as geometric
from check_raar_coverings_ad import Dual, Model, exact
from search_raar_circles import step


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    mp.iv.dps = 100
    mp.mp.dps = 120
    root = Path(__file__).resolve().parents[1]
    source = root/'results/raar-circles/render/case-1051-animation.npz'
    data = np.load(source)
    states, p = data['states'], data['params']
    trig = np.array([[np.nextafter(float(v.a), -np.inf), np.nextafter(float(v.b), np.inf)]
                     for v in [mp.iv.sin(exact(p[2])), mp.iv.cos(exact(p[2]))]])
    for interval, v in zip(trig, [mp.iv.sin(exact(p[2])), mp.iv.cos(exact(p[2]))]):
        assert exact(interval[0]) < v.a and exact(interval[1]) > v.b
    rng = np.random.default_rng(202611011033)
    primitive_checks = 0
    for _ in range(2000):
        values = np.ldexp(rng.uniform(-1, 1, 4), rng.integers(-100, 101, 4))
        a, b = np.sort(values[:2]), np.sort(values[2:])
        fa, fb = [Fraction(float(v)) for v in a], [Fraction(float(v)) for v in b]
        for operation, answers in [(fast.add, [fa[0]+fb[0], fa[1]+fb[1]]),
                                   (fast.mul, [x*y for x in fa for y in fb])]:
            result = operation(a, b)
            assert Fraction(float(result[0])) <= min(answers)
            assert Fraction(float(result[1])) >= max(answers)
            primitive_checks += 1
        if not a[0] <= 0 <= a[1]:
            result = fast.reciprocal(a)
            assert Fraction(float(result[0])) <= 1/fa[1]
            assert Fraction(float(result[1])) >= 1/fa[0]
            primitive_checks += 1
        a = np.sort(np.abs(a))
        result = fast.root(a)
        assert mp.mpf(float(result[0])) <= mp.sqrt(mp.mpf(float(a[0])))
        assert mp.mpf(float(result[1])) >= mp.sqrt(mp.mpf(float(a[1])))
        primitive_checks += 1
    print(json.dumps({'primitive_checks': primitive_checks, 'stage': 'compiling map'}), flush=True)
    model = Model(p)
    rows, point_checks = [], 0
    enclosure = geometric.remainder if args.method == 'geometric' else fast.centered_remainder
    sample_ids = np.linspace(10000, len(states)-1, 64, dtype=int)
    for radius in [0., .0001, .00125, .005]:
        widths, failures = [], 0
        for ident in sample_ids:
            center = states[ident]
            box = np.column_stack((np.nextafter(center-radius, -np.inf),
                                   np.nextafter(center+radius, np.inf)))
            values, jac, axis, flag = fast.evaluate_box(box, p, trig)
            if flag:
                failures += 1
                continue
            image, linear, _, numeric_flag = step(center, p)
            assert numeric_flag == 0
            radii, remainder, _, flag = enclosure(box, center, image, linear, p, trig, 2.**-35)
            assert flag == 0 and np.isfinite(remainder).all()
            for point in [center, box[:, 0], box[:, 1], rng.uniform(box[:, 0], box[:, 1])]:
                x = [Dual(exact(point[i]), [mp.iv.mpf(int(i == j)) for j in range(3)]) for i in range(3)]
                actual, actual_axis = model.field(x)
                assert exact(axis) < actual_axis
                for i in range(3):
                    assert exact(values[i, 0]) <= actual[i].v.a and exact(values[i, 1]) >= actual[i].v.b
                    for j in range(3):
                        assert exact(jac[i, j, 0]) <= actual[i].d[j].a and exact(jac[i, j, 1]) >= actual[i].d[j].b
                    affine = exact(image[i]) + sum(exact(linear[i, j])*(exact(point[j])-exact(center[j])) for j in range(3))
                    residual = actual[i].v-affine
                    assert residual.a > -exact(remainder[i]) and residual.b < exact(remainder[i])
                point_checks += 1
            widths.append(remainder.tolist())
        rows.append({'radius': radius, 'boxes': len(sample_ids), 'extension_failures': failures,
                     'maximum_remainder': np.max(widths, axis=0).tolist() if widths else None})
        print(json.dumps(rows[-1]), flush=True)
    begin = time.time()
    for ident in sample_ids:
        center = states[ident]
        fast.evaluate_box(np.column_stack((center-.00125, center+.00125)), p, trig)
    begin_enclosure = time.time()
    for ident in sample_ids:
        center = states[ident]
        image, linear, _, _ = step(center, p)
        enclosure(np.column_stack((center-.00125, center+.00125)), center, image, linear, p, trig, 2.**-35)
    report = {'method': args.method, 'primitive_checks': primitive_checks, 'independent_point_and_jacobian_checks': point_checks,
              'rows': rows, 'parameters': p.tolist(), 'trig_intervals': trig.tolist(),
              'seconds_per_enclosure': (time.time()-begin_enclosure)/len(sample_ids), 'total_seconds': time.time()-started,
              'reference_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'scope': 'Implementation sanity checks, not exhaustive tests or a trapping certificate. Arithmetic enclosures rely on the stated IEEE binary64 and correctly rounded sqrt assumptions.'}
    for name in ['raar_interval64.py', 'raar_geometric_enclosure.py', 'validate_raar_interval64.py']:
        path = Path(__file__).with_name(name)
        (args.out/name).write_bytes(path.read_bytes())
        report[name+'_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--method', choices=['derivative', 'geometric'], default='derivative')
    main(parser.parse_args())

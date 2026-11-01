"""Exact rational audit of orientation partitions and periodic certificates."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import json
from pathlib import Path


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def endpoint(t):
    sign, mantissa, exponent, bits = t
    assert sign in [0, 1] and bits == mantissa.bit_length()
    return (-1 if sign else 1)*Q(mantissa)*Q(2)**exponent


def interval(row):
    result = tuple(map(endpoint, row['binary']))
    assert result[0] <= result[1]
    return result


def point(x):
    return Q(x), Q(x)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    v = [x*y for x in a for y in b]
    return min(v), max(v)


def isum(v):
    total = point(0)
    for a in v:
        total = add(total, a)
    return total


def contains(a, b):
    return a[0] <= b[0] <= b[1] <= a[1]


def nonzero(a):
    return a[0] > 0 or a[1] < 0


def opposite(a, b):
    return a[1] < 0 < b[0] or b[1] < 0 < a[0]


def mm(a, b):
    return [[isum(mul(a[i][k], b[k][j]) for k in range(3)) for j in range(3)] for i in range(3)]


def mv(a, x):
    return [isum(mul(a[i][j], x[j]) for j in range(3)) for i in range(3)]


def det(a):
    return add(sub(mul(a[0][0], sub(mul(a[1][1], a[2][2]), mul(a[1][2], a[2][1]))),
                   mul(a[0][1], sub(mul(a[1][0], a[2][2]), mul(a[1][2], a[2][0])))),
               mul(a[0][2], sub(mul(a[1][0], a[2][1]), mul(a[1][1], a[2][0]))))


def source_check(folder):
    report = json.loads((folder/'report.json').read_text(encoding='utf-8'))
    assert report['status'] == 'passed'
    for name, digest in report['sources'].items():
        assert sha(folder/name) == digest
    return report


def orientation(folder, windows):
    report = source_check(folder)
    assert report['input_sha256'] == sha(windows) == sha(folder/'input.json')
    total, rows = point(1), []
    assert [r['window'] for r in report['windows']] == list(range(23))
    for r in report['windows']:
        path = folder/r['file']
        assert sha(path) == r['sha256']
        tree, count, volume, bounds, min_radius = {}, 0, Q(0), [], None
        with gzip.open(path, 'rt', encoding='utf-8') as handle:
            for line in handle:
                row = json.loads(line)
                address, node = row['path'], tree
                assert len(address) % 2 == 0
                for k in range(0, len(address), 2):
                    token = address[k:k+2]
                    assert token[0] in '012' and token[1] in '01' and '#' not in node
                    node = node.setdefault(token, {})
                assert not node
                node['#'] = True
                count += 1
                volume += Q(1, 2**(len(address)//2))
                value, product = interval(row['determinant']), point(1)
                assert len(row['one_step']) == len(row['radii']) == 4
                for d in row['one_step']:
                    product = mul(product, interval(d))
                assert contains(value, product)
                assert value[1] < 0 if r['sign'] == -1 else value[0] > 0
                for radius in row['radii']:
                    rr = interval(radius)[0]
                    assert rr > 0
                    min_radius = rr if min_radius is None else min(min_radius, rr)
                bounds.append(value)
        def visit(node):
            if '#' in node:
                assert set(node) == {'#'}
                return 1
            keys = sorted(node)
            assert len(keys) == 2 and keys[0][0] == keys[1][0] and keys[0][1] == '0' and keys[1][1] == '1'
            return sum(visit(v) for v in node.values())
        assert visit(tree) == count == r['leaves'] and volume == 1
        hull = min(a for a, b in bounds), max(b for a, b in bounds)
        assert contains(interval(r['determinant']), hull)
        total = mul(total, interval(r['determinant']))
        rows.append(dict(window=r['window'], cells=count, sign=r['sign'], minimum_radius=float(min_radius)))
    assert contains(interval(report['excursion_determinant']), total) and total[1] < 0
    assert [r['window'] for r in rows if r['sign'] == -1] == [4, 5, 6]
    proper_divisors = [d for d in range(1, 92) if 92 % d == 0]
    assert proper_divisors == [1, 2, 4, 23, 46]
    assert all(4 % d == 0 for d in [1, 2, 4])
    assert all((92//d) % 2 == 0 for d in [23, 46])
    return dict(source_sha256=sha(folder/'report.json'), cells=sum(r['cells'] for r in rows), windows=rows,
                determinant_outer=[str(total[0]), str(total[1])], determinant_display=[float(total[0]), float(total[1])],
                least_period=92, divisor_exclusions={'N0_N1_separation': [1, 2, 4], 'negative_determinant_even_power': [23, 46]})


def cycles(folder, windows):
    report = source_check(folder)
    assert report['windows_sha256'] == sha(windows) == sha(folder/'windows.json')
    raw = json.loads(windows.read_text(encoding='utf-8'))
    rows = []
    I = [[point(int(i == j)) for j in range(3)] for i in range(3)]
    for row in report['cycles']:
        box = list(map(interval, row['box']))
        center = list(map(point, row['center_decimal']))
        C = [[point(v) for v in r] for r in row['preconditioner_decimal']]
        assert nonzero(det(C))
        M = [[interval(v) for v in r] for r in row['monodromy']]
        H = [[sub(M[i][j], I[i][j]) for j in range(3)] for i in range(3)]
        CH = mm(C, H)
        R = [[sub(I[i][j], CH[i][j]) for j in range(3)] for i in range(3)]
        contraction = max(sum(max(abs(v[0]), abs(v[1])) for v in r) for r in R)
        assert contraction < 1
        residual = list(map(interval, row['residual_at_center']))
        correction = mv(C, residual)
        remainder = mv(R, [sub(a, b) for a, b in zip(box, center)])
        K = [add(sub(a, b), c) for a, b, c in zip(center, correction, remainder)]
        assert all(b[0] < k[0] <= k[1] < b[1] for b, k in zip(box, K))
        assert all(contains(interval(r), k) for r, k in zip(row['krawczyk'], K))
        assert contains(interval(row['determinant']), det(M))
        period = row['period']
        divisors = [d for d in range(1, period) if period % d == 0]
        assert [r['divisor'] for r in row['proper_divisor_exclusions']] == divisors
        for d in divisors:
            assert any(nonzero(sub(interval(v), b)) for v, b in zip(row['orbit'][d], box))
        assert all(interval(v)[0] > 0 for v in row['radii'])
        coefficients = [neg(isum(M[i][i] for i in range(3))),
                        isum(sub(mul(M[i][i], M[j][j]), mul(M[i][j], M[j][i])) for i in range(3) for j in range(i+1, 3)), neg(det(M))]
        def polynomial(x):
            return add(mul(add(mul(add(x, coefficients[0]), x), coefficients[1]), x), coefficients[2])
        brackets = []
        for mult in row['multipliers']:
            lo, hi = map(Q, mult['endpoints_decimal'])
            assert lo < hi and opposite(polynomial(point(lo)), polynomial(point(hi)))
            assert mult['stable'] == (-1 < lo < hi < 1)
            assert mult['unstable'] == (hi < -1 or lo > 1)
            brackets.append([lo, hi])
        assert len(brackets) == 3 and all(brackets[i][1] < brackets[i+1][0] for i in range(2))
        local = mv([[point(v) for v in r] for r in raw['exact_windows'][0]['inverse']], [sub(box[j], point(raw['centers'][0][j])) for j in range(3)])
        membership = all(-1 < a <= b < 1 for a, b in local)
        assert membership == row['in_window_zero']
        rows.append(dict(id=row['id'], period=period, contraction=float(contraction), in_window_zero=membership,
                         multipliers=[[float(v) for v in b] for b in brackets], determinant=[float(v) for v in det(M)]))
    return dict(source_sha256=sha(folder/'report.json'), cycles=rows)


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    coverage = json.loads(args.covering_audit.read_text(encoding='utf-8'))
    assert coverage['status'] == 'passed' and coverage['complete_graph'] and coverage['windows'] == 23
    assert coverage['input_sha256'] == sha(args.windows)
    separation = list(map(Q, coverage['exact_separation']))
    assert separation[1] < -1 or separation[0] > 1
    result = dict(status='passed', windows_sha256=sha(args.windows), covering_audit_sha256=sha(args.covering_audit),
                  orientation=[orientation(p, args.windows) for p in args.orientation],
                  cycles=[cycles(p, args.windows) for p in args.cycles],
                  scope='Exact binary endpoints, complete partitions, determinant products, least-period logic, rational Krawczyk arithmetic and real-root brackets. Function and derivative enclosures are supplied by two formula implementations sharing mpmath.')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status='passed', orientation=[r['cells'] for r in result['orientation']], cycles=[len(r['cycles']) for r in result['cycles']])))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['windows', 'covering-audit', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--orientation', type=Path, nargs='+', required=True)
    p.add_argument('--cycles', type=Path, nargs='+', required=True)
    main(p.parse_args())

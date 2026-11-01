"""Exact binary-endpoint and complete-tree audit for product-map coverings."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import gzip
import hashlib
import json
from pathlib import Path


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def binary(t):
    sign, mantissa, exponent, bits = t
    assert sign in [0, 1] and mantissa >= 0 and bits == mantissa.bit_length()
    return (-1 if sign else 1)*Q(mantissa)*Q(2)**exponent


def partition(path, kind, orientation):
    tree, depths, modes = {}, Counter(), Counter()
    minimum, radius, count = None, None, 0
    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            row = json.loads(line)
            address, node = row['path'], tree
            assert len(address) % 2 == 0
            for k in range(0, len(address), 2):
                token = address[k:k+2]
                assert token[0] in '012' and token[1] in '01'
                assert kind == 'volume' or token[0] != '0'
                assert '#' not in node
                node = node.setdefault(token, {})
            assert not node
            node['#'] = True
            count += 1
            depths[len(address)//2] += 1
            image = [[binary(t) for t in v] for v in row['image_binary']]
            assert all(a <= b for a, b in image)
            rr = binary(row['radius_binary'])
            assert rr > 0
            radius = rr if radius is None else min(radius, rr)
            if kind == 'volume':
                if image[0][0] > 1:
                    margin, mode = image[0][0]-1, 'right-exterior'
                elif image[0][1] < -1:
                    margin, mode = -1-image[0][1], 'left-exterior'
                else:
                    margin = min(v for a, b in image[1:] for v in [a+1, 1-b])
                    mode = 'stable-interior'
            else:
                sign = -orientation if kind == 'left' else orientation
                margin = (image[0][0] if sign == 1 else -image[0][1])-1
                mode = 'exit-face'
            assert margin > 0
            minimum = margin if minimum is None else min(minimum, margin)
            modes[mode] += 1
    def visit(node):
        if '#' in node:
            assert set(node) == {'#'}
            return 1
        keys = sorted(node)
        assert len(keys) == 2 and keys[0][0] == keys[1][0] and keys[0][1] == '0' and keys[1][1] == '1'
        return sum(visit(v) for v in node.values())
    assert visit(tree) == count and sum(Q(n, 2**d) for d, n in depths.items()) == 1
    return dict(cells=count, maximum_depth=max(depths), minimum_radius_exact=str(radius), minimum_radius=float(radius), minimum_margin_exact=str(minimum), minimum_margin=float(minimum), modes=dict(modes))


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    report = json.loads((args.source/'report.json').read_text())
    raw = json.loads((args.source/'input.json').read_text())
    assert report['status'] == 'passed' and sha(args.source/'input.json') == report['input_sha256']
    for name, digest in report['sources'].items():
        assert sha(args.source/name) == digest
    size = len(raw['centers'])
    assert raw['edges'] == [[0, 0], [0, 1]]+[[i, (i+1) % size] for i in range(1, size)]
    assert len(raw['orientation']) == size+1
    matrices, inverses = [], []
    for b, r, stored in zip(raw['bases'], raw['radii'], raw['exact_windows']):
        s = [[Q(v)*Q(r[j]) for j, v in enumerate(row)] for row in b]
        t = [[Q(v) for v in row] for row in stored['inverse']]
        assert s == [[Q(v) for v in row] for row in stored['matrix']]
        assert all(sum(s[i][k]*t[k][j] for k in range(3)) == int(i == j) for i in range(3) for j in range(3))
        matrices.append(s)
        inverses.append(t)
    d = [Q(b)-Q(a) for a, b in zip(raw['centers'][0], raw['centers'][1])]
    center = sum(inverses[0][0][k]*d[k] for k in range(3))
    radius = sum(abs(sum(inverses[0][0][k]*matrices[1][k][j] for k in range(3))) for j in range(3))
    assert center+radius < -1 or center-radius > 1
    checks, seen = [], set()
    for check in report['checks']:
        e, kind = check['edge'], check['kind']
        assert e in range(size+1) and kind in ['volume', 'left', 'right'] and (e, kind) not in seen
        seen.add((e, kind))
        assert check['status'] == 'passed' and check['pending'] == 0
        path = args.source/f'edge-{e}-{kind}.jsonl.gz'
        assert sha(path) == check['certificate_sha256']
        row = partition(path, kind, raw['orientation'][e])
        assert row['cells'] == check.get('leaves', check.get('accepted'))
        assert row['maximum_depth'] == check['maximum_depth']
        assert 2*row['cells']-1 == check['examined']
        checks.append(dict(edge=e, kind=kind, **row))
    assert len(seen) == 3*(size+1)
    log_lower = 2*sum(Q(1, (2*k+1)*3**(2*k+1)) for k in range(12))
    result = dict(status='passed', source_report_sha256=sha(args.source/'report.json'), input_sha256=report['input_sha256'], windows=size, complete_graph=True, total_leaves=sum(r['cells'] for r in checks), checks=checks, exact_separation=[str(center-radius), str(center+radius)], minimum_margin=min(r['minimum_margin'] for r in checks), minimum_radius=min(r['minimum_radius'] for r in checks), entropy_bound=f'log(2)/{4*size} per F iterate on a compact invariant subsystem', entropy_rational_lower_bound=str(log_lower/(4*size)), entropy_decimal_lower_bound=float(log_lower/(4*size)), scope='Exact binary endpoint, complete tiling, strict inequality, matrix inverse, separation and hash audit. Enclosures need independent recomputation. No trapping or attraction claim.')
    (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['checks', 'exact_separation']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['source', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    main(p.parse_args())

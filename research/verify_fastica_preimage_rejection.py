"""Audit the retained preimage certificate and its complete box partition."""
from pathlib import Path
from fractions import Fraction
import hashlib
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[1]
PACKET = ROOT/'results/fastica-preimage-rejection'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_partition(leaves):
    """Reconstruct both complete binary trees from leaf addresses alone."""
    counts = {'splits': 0, 'leaves': 0, 'maximum_depth': 0}
    seen = {(row['chart'], row['address']) for row in leaves}
    assert len(seen) == len(leaves)
    assert all(chart in [0, 1] and len(code) % 2 == 0 and
               all(code[i] in '012' and code[i+1] in '01' for i in range(0, len(code), 2))
               for chart, code in seen)

    def visit(chart, prefix):
        descendants = [code for c, code in seen if c == chart and code.startswith(prefix)]
        assert descendants, (chart, prefix)
        if prefix in descendants:
            assert len(descendants) == 1, 'Overlapping leaf and descendant'
            counts['leaves'] += 1
            counts['maximum_depth'] = max(counts['maximum_depth'], len(prefix)//2)
            return
        children = {code[len(prefix):len(prefix)+2] for code in descendants}
        coordinates = {code[0] for code in children}
        assert len(coordinates) == 1
        coordinate = next(iter(coordinates))
        assert children == {coordinate+'0', coordinate+'1'}, 'Missing child'
        counts['splits'] += 1
        for bit in '01':
            visit(chart, prefix+coordinate+bit)

    for chart in [0, 1]:
        visit(chart, '')
    for row in leaves:
        if row['reason'] == 'outside_tube':
            assert Fraction(row['positive_gap']) > 0
        else:
            assert row['reason'] == 'excluded_equation' and row['equation'] in [0, 1, 2]
            lo, hi = map(Fraction, row['interval'])
            assert lo <= hi and (hi < 0 or lo > 0)
    assert counts['leaves'] == counts['splits'] + 2
    return counts


def main():
    manifest = json.loads((PACKET/'manifest.json').read_text())
    for row in manifest['files']:
        path = PACKET/row['path']
        assert path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], row['path']
    for row in manifest['dependencies']:
        assert sha(ROOT/row['path']) == row['sha256'], row['path']
    report = json.loads((PACKET/'report.json').read_text())
    spec = importlib.util.spec_from_file_location('preimage_certificate', ROOT/'research/fastica_preimage_rejection.py')
    check = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(check)
    target = [check.Interval(c-check.RADIUS, c+check.RADIUS) for c in check.CENTER]
    assert report['status'] == 'passed' and Fraction(report['epsilon']) == check.E
    assert Fraction(report['tube_squared_ratio']) == check.K2
    assert list(map(Fraction, report['previous_raw_vector'])) == check.PREVIOUS
    assert report['target_box'] == [v.pair() for v in target]
    assert report['source_sha256'] == sha(ROOT/'research/fastica_preimage_rejection.py') == sha(PACKET/'check-source.py')
    assert report['root_certificate'] == check.certify_root(target)
    assert report['predecessor_exclusion'] == check.exclude_predecessor(target)
    root = report['root_certificate']
    assert Fraction(root['contraction_infinity_bound']) < Fraction('0.000226')
    assert min(map(Fraction, root['strict_self_map_margins'])) > Fraction('9.70e-8')
    assert Fraction(root['previous_tube_margin']) > Fraction('0.1178')
    assert Fraction(root['root_box_tube_margin']) > Fraction('0.1055')
    assert Fraction(root['separation_in_C_over_minus_B']) > Fraction('0.0791')
    accounting = audit_partition(report['predecessor_exclusion']['leaves'])
    assert accounting == {'splits': 43, 'leaves': 45, 'maximum_depth': 7}
    exploration = json.loads((PACKET/'exploration.json').read_text())
    assert exploration['source_sha256'] == sha(PACKET/'exploration-source.py')
    assert exploration['model_sha256'] == sha(ROOT/'results/fastica-tripling/code/model.py')
    assert len(exploration['rows']) == 16 and sum(len(r['attempts']) for r in exploration['rows']) == 448
    assert [len(r['roots']) for r in exploration['rows']] == [2 if i == 13 else 1 for i in range(16)]
    assert list(map(lambda v: Fraction(str(v)), exploration['rows'][13]['known_previous'])) == check.PREVIOUS
    print(json.dumps({'status': 'passed', 'packet_files': len(manifest['files']),
                      'partition': accounting, 'exploratory_solves': 448,
                      'scope': 'Independent address accounting; arithmetic replay shares the certificate implementation.'}))


if __name__ == '__main__':
    main()

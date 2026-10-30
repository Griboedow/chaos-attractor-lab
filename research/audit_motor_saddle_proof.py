"""Audit retained saddle certificates and provenance without reintegration."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output; preserve prior evidence')
    root = Path(__file__).resolve().parents[1]
    folder = root / 'results/motor-saddle-proof'
    source_files = list(folder.rglob('*.py')) + list((root / 'research').glob('*.py'))
    source_index = {}
    for path in source_files:
        source_index.setdefault(digest(path), []).append(path.relative_to(root).as_posix())
    checks = []

    def check(name, passed, **details):
        checks.append(dict(name=name, passed=bool(passed), **details))

    def read(name):
        return json.loads((folder / name).read_text(encoding='utf-8'))

    for path in sorted(folder.rglob('*.json')):
        data = json.loads(path.read_text(encoding='utf-8'))
        hashes = data.get('source_sha256', {})
        if isinstance(hashes, dict):
            for name, expected in hashes.items():
                matches = source_index.get(expected, [])
                check('retained executed source', bool(matches), record=path.relative_to(root).as_posix(), source=name, sha256=expected, matches=matches)
        if 'parent_certificate' in data:
            parent = data['parent_certificate']
            check('postprocessor parent hash', digest(root / parent['path']) == parent['sha256'], record=path.name)
        if 'postprocessor_source_snapshot' in data:
            current = path.parent / data['postprocessor_source_snapshot']
            check('executing postprocessor snapshot', digest(current) == data['postprocessor_source_sha256'], record=path.name)
        for key in ('center_input', 'input_path'):
            expected_key = 'center_input_sha256' if key == 'center_input' else 'input_sha256'
            if key in data and expected_key in data:
                check('input hash', digest(root / data[key]) == data[expected_key], record=path.name, input=data[key])

    for name in ('similarity-certificate.json', 'recheck-01.json', 'recheck-postprocess-01.json'):
        data = read(name)
        check('full saddle status', data['status'] == 'six_cluster_saddle_certified', record=name)
        check('existence gate', data['krawczyk_passed'], record=name)
        check('six distinct clusters', len(data['pairwise_cluster_distinctions']) == 15 and all(x['distinct'] for x in data['pairwise_cluster_distinctions']), record=name)
        check('complete Floquet count', data['spectrum']['floquet_count'] == dict(unstable=1, neutral=1, stable=21), record=name)
        check('similarity Schur gate', data['spectrum']['similarity_certificate']['passed'], record=name)
        basis = data['spectrum']['similarity_certificate']
        check('exact basis file hash', digest(root / basis['basis_path']) == basis['basis_sha256'], record=name)
        check('both contrast gates', all(x['test']['passed'] for x in data['spectrum']['transverse'].values()), record=name)

    period = read('least-period.json')
    check('least-period status', period['status'] == 'least_saddle_and_row_periods_certified')
    for name, value in period['gates'].items():
        check('least-period gate: ' + name, value)
    check('four upward and downward full-state crossings', sum(x['upward'] for x in period['full_crossings']) == 4 and sum(x['downward'] for x in period['full_crossings']) == 4)
    check('two upward and downward row crossings', period['row2_crossings']['upward'] == 2 and period['row2_crossings']['downward'] == 2)
    post = read('recheck-postprocess-01.json')
    old_hash = post['source_sha256']['certify_motor_saddle_spectrum.py']
    check('different helper versions preserved', old_hash != post['postprocessor_source_sha256'] and digest(folder / 'recheck-postprocess-01.certify_motor_saddle_spectrum.py') == old_hash)
    for name in ('analytic-controls-01.json', 'independent-review-controls.json', 'independent-similarity-controls/analysis.json', 'zero-window-controls.json'):
        check('implementation control status', read(name)['status'] == 'passed', record=name)
    check('inconclusive raw quotient result retained', read('attempt-02.json')['status'] == 'periodic_six_cluster_orbit_certified')
    files = sorted(p for p in folder.rglob('*') if p.is_file())
    files += [root / 'research' / name for name in ('validate_motor_saddle_cycle.py', 'validate_motor_saddle_period.py', 'certify_motor_saddle_spectrum.py', 'check_motor_saddle_proof.py', 'audit_motor_saddle_proof.py')]
    files += [root / 'docs/motor-saddle-certificate.md']
    result = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='passed' if all(x['passed'] for x in checks) else 'failed', checks=checks,
        files=[dict(path=p.relative_to(root).as_posix(), bytes=p.stat().st_size, sha256=digest(p)) for p in files],
        limitations=['Stored-gate and provenance audit only; no new integration or independent proof implementation.'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(result['status'], len(checks), 'checks', len(files), 'files')
    if result['status'] != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()

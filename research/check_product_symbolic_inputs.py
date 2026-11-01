"""Bind the symbolic certificate to the exact published Case 988 inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    raw = json.loads((args.certificate/'input.json').read_text(encoding='utf-8'))
    first = json.loads((args.certificate/'report.json').read_text(encoding='utf-8'))
    second = json.loads((args.independent/'report.json').read_text(encoding='utf-8'))
    assert first['status'] == second['status'] == 'passed'
    assert first['input_sha256'] == second['input_sha256'] == sha(args.certificate/'input.json') == sha(args.independent/'input.json')
    assert second['primary_report_sha256'] == sha(args.certificate/'report.json')
    for folder, source in [(args.primary_audit, args.certificate), (args.independent_audit, args.independent)]:
        audit = json.loads((folder/'report.json').read_text(encoding='utf-8'))
        assert audit['status'] == 'passed' and audit['source_report_sha256'] == sha(source/'report.json')
    with np.load(args.dense) as data:
        assert data['params'].tolist() == raw['params']
        assert data['points'].tolist() == raw['points']
        assert int(data['order']) == 2
        states = data['states']
    distances, indices = cKDTree(states).query(np.array(raw['centers']))
    files = [args.dense, args.certificate/'input.json', args.certificate/'report.json', args.independent/'report.json', args.primary_audit/'report.json', args.independent_audit/'report.json']
    result = dict(status='passed', case=988, dependencies={str(p): sha(p) for p in files}, exact_parameters_hex=[float(v).hex() for v in raw['params']], exact_sites_hex=[[float(v).hex() for v in p] for p in raw['points']], same_exact_instance=True, window_centers_to_finite_cloud=dict(states=len(states), distances=distances.tolist(), nearest_state_indices=indices.tolist(), maximum=float(distances.max()), scope='Numerical distances to a finite cloud, not membership in its infinite closure or an attracting set.'), scope='Bitwise model-input and certificate-dependency binding. Covering enclosures and complete trees are verified by their separate programs.')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status='passed', same_exact_instance=True, maximum_center_distance=float(distances.max()))))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['dense', 'certificate', 'independent', 'primary-audit', 'independent-audit', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    main(p.parse_args())

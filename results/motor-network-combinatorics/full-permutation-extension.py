"""Enumerate all labeled saddle images and the formal S8 departure graph.

Uses integer permutations and exact rational checks on archived certificate
gates. It does not repeat interval integration or prove connecting orbits.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
from itertools import combinations, permutations
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
I = tuple(range(8))


def compose(a, b):
    return tuple(a[b[j]] for j in range(8))


def inverse(a):
    return tuple(a.index(j) for j in range(8))


def swap(a, b):
    result = list(I)
    result[a], result[b] = result[b], result[a]
    return tuple(result)


def closure(generators):
    group, pending = {I}, [I]
    while pending:
        a = pending.pop()
        for b in generators:
            p = compose(a, b)
            if p not in group:
                group.add(p)
                pending.append(p)
    assert all(inverse(a) in group for a in group)
    return group


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dyadic(pair):
    n, e = int(pair[0]), int(pair[1])
    return Fraction(n * 2**e) if e >= 0 else Fraction(n, 2**(-e))


def bounds(record):
    m, r = dyadic(record['midpoint']), dyadic(record['radius'])
    return m-r, m+r


def strict_schur(test):
    return (test['passed'] and bounds(test['radius'])[1] < 1
            and all(bounds(step['margin'])[0] > 0 for step in test['steps']))


def certificate_checks(period_path, spectral_path):
    period = json.loads(period_path.read_text())
    spectral = json.loads(spectral_path.read_text())
    checks = {}
    checks['period_status'] = period['status'] == 'least_saddle_and_row_periods_certified'
    checks['all_period_gates'] = all(period['gates'].values())
    checks['period_links_same_spectral_certificate'] = period['input_sha256'] == sha(spectral_path)
    checks['anchor_links_same_certificate'] = period['anchor_revalidation']['center_input_sha256'] == sha(spectral_path)
    checks['period_and_row_period_positive'] = min(bounds(period['positive_period'])[0],
        bounds(period['row2_least_period'])[0]) > 0
    checks['complete_saddle_status'] = spectral['status'] == 'six_cluster_saddle_certified'
    checks['nonconstant_orbit'] = spectral['krawczyk_passed'] and bounds(spectral['section_speed'])[0] > 0
    distinctions = spectral['pairwise_cluster_distinctions']
    checks['all_fifteen_cluster_distinctions'] = (
        {tuple(v['rows']) for v in distinctions} == set(combinations(range(6), 2))
        and all(v['distinct'] and any(bounds(x)[0] > 0 or bounds(x)[1] < 0
            for x in v['differences']) for v in distinctions))
    checks['exact_cluster_layout_and_half_symmetry'] = (
        spectral['multiplicities'] == [1,1,2,1,1,2]
        and spectral['half_permutation'] == [1,0,2,4,3,5])
    stable = spectral['spectrum']['transverse']['2']['test']
    expanding = spectral['spectrum']['transverse']['5']['test']
    checks['pair_23_strictly_stable'] = strict_schur(stable)
    checks['pair_67_has_expanding_real_root'] = (
        expanding['passed'] and Fraction(expanding['expanding_root_bracket'][0]) > 1
        and bounds(expanding['left'])[1] < 0 and bounds(expanding['right'])[0] > 0)
    checks['other_pair_67_roots_stable'] = strict_schur(expanding['remaining_quadratic_schur'])
    snapshots = {}
    for path, record in ((period_path, period), (spectral_path, spectral)):
        for name, digest in record['source_sha256'].items():
            snapshot = path.with_name(path.stem+'.'+name)
            snapshots[snapshot.relative_to(ROOT).as_posix()] = sha(snapshot) == digest
    postprocessor = spectral_path.with_name(spectral_path.stem+'.certify_motor_saddle_spectrum.py')
    snapshots[postprocessor.relative_to(ROOT).as_posix()] = sha(postprocessor) == spectral['postprocessor_source_sha256']
    checks['all_certificate_source_snapshots_match'] = all(snapshots.values())
    assert all(checks.values()), checks
    return dict(passed=True, checks=checks, source_snapshot_checks=snapshots,
        least_period=period['positive_period'], row2_least_period=period['row2_least_period'],
        exact_rho_decimal=spectral['exact_rho_decimal'],
        pair_23_schur_radius=stable['radius'], pair_67_expanding_bracket=expanding['expanding_root_bracket'],
        limitation='Checks archived theorem gates and source provenance; no interval integration or monodromy computation is repeated.')


def enumerate_graph(reference):
    h0 = closure((swap(2,3), swap(6,7)))
    half = (1,0,3,2,5,4,7,6)
    h = closure((*h0, half))
    plus, minus = (tuple(reference['generators'][name]) for name in ('plus','minus'))
    g = closure((*h, plus, minus))
    assert h == {tuple(p) for p in reference['known_orbit_subgroup']}
    assert g == {tuple(p) for p in reference['group']}
    assert len(h0) == 4 and len(h) == 8 and len(g) == 64
    all_permutations = list(permutations(range(8)))
    node_of, nodes = {}, []
    for a in all_permutations:
        if a in node_of:
            continue
        coset = {compose(b,a) for b in h}
        assert len(coset) == 8 and min(coset) == a
        index = len(nodes)
        nodes.append(dict(id=index, representative=a))
        for b in coset:
            assert b not in node_of
            node_of[b] = index
    component_of, components, local_name = {}, [], {}
    for a in all_permutations:
        if a in component_of:
            continue
        coset = {compose(b,a) for b in g}
        assert len(coset) == 64 and min(coset) == a
        index = len(components)
        mapping = {item['id']:node_of[compose(tuple(item['representative']),a)]
                   for item in reference['nodes']}
        assert len(set(mapping.values())) == 8
        assert set(mapping.values()) == {node_of[b] for b in coset}
        components.append(dict(id=index, representative=a, nodes_by_reference_label=mapping))
        for name, node in mapping.items():
            local_name[node] = name
        for b in coset:
            assert b not in component_of
            component_of[b] = index
    representative_checks = 0
    for node in nodes:
        a = node['representative']
        destinations = {node_of[compose(b,a)] for b in (plus,minus)}
        assert len(destinations) == 2
        for b in h:
            representative_checks += 1
            assert {node_of[compose(departure,compose(b,a))]
                    for departure in (plus,minus)} == destinations
        node.update(component=component_of[a], local_label=local_name[node['id']],
                    outgoing=sorted(destinations))
    indegrees = Counter(target for node in nodes for target in node['outgoing'])
    assert set(indegrees.values()) == {2}
    for component in components:
        members = set(component['nodes_by_reference_label'].values())
        for start in members:
            pending, reached = [start], set()
            while pending:
                node = pending.pop()
                if node in reached:
                    continue
                reached.add(node)
                assert nodes[node]['component'] == component['id']
                pending += [target for target in nodes[node]['outgoing'] if target not in reached]
            assert reached == members
        # Also verify the full eight-node adjacency agrees with the reference
        # under right translation by the component representative.
        mapping = component['nodes_by_reference_label']
        for local in reference['nodes']:
            assert set(nodes[mapping[local['id']]]['outgoing']) == {mapping[name] for name in local['outgoing']}
    pair_sets = {frozenset((2,3)), frozenset((6,7))}
    normalizer = [a for a in all_permutations if
        {frozenset(a[j] for j in pair) for pair in pair_sets} == pair_sets]
    pair_preserving = [a for a in normalizer if all(frozenset(a[j] for j in pair) == pair for pair in pair_sets)]
    assert all({compose(compose(a,b),inverse(a)) for b in h0} == h0 for a in normalizer)
    assert len(normalizer) == 192 and len(pair_preserving) == 96
    assert len(all_permutations) == len(node_of) == len(component_of) == 40320
    assert len(nodes) == 5040 and len(components) == 630
    return dict(group_order=40320, pointwise_isotropy_order=4, orbital_isotropy_order=8,
        pointwise_isotropy=sorted(h0), orbital_isotropy=sorted(h),
        pointwise_normalizer_order=192, pair_preserving_normalizer_order=96,
        departure_subgroup_order=64, node_count=len(nodes), edge_count=2*len(nodes),
        strongly_connected_component_count=len(components), nodes_per_component=8,
        edges_per_component=16, all_node_indegrees=2, all_node_outdegrees=2,
        coset_representative_checks=representative_checks, strong_connectivity_start_checks=len(nodes),
        action='(q dot a)[j]=q[a[j]]; compose(a,b)[j]=a[b[j]]',
        nodes='Left cosets H a; each listed representative is lexicographically least.',
        edge_rule='H a -> H compose(g_plus,a), H compose(g_minus,a).',
        components_convention='Left cosets G a. Right translation by a identifies each component with the reference eight-node graph.',
        node_records=nodes, component_records=components)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve prior evidence: choose a new output')
    started = time.monotonic()
    period = ROOT/'results/motor-saddle-proof/least-period.json'
    spectral = ROOT/'results/motor-saddle-proof/similarity-certificate.json'
    reference = ROOT/'results/motor-network-combinatorics/exact-graph.json'
    report = ROOT/'docs/motor-saddle-isotropy.md'
    gates = certificate_checks(period, spectral)
    graph = enumerate_graph(json.loads(reference.read_text()))
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        status='exact_full_permutation_extension_passed', wall_seconds=time.monotonic()-started,
        source_sha256=sha(Path(__file__)), report_sha256=sha(report),
        input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in (period,spectral,reference)},
        certificate_checks=gates, full_graph=graph,
        conclusion='The linked interval theorem inputs and the analytic isotropy argument give exactly 5040 distinct labeled periodic-orbit images. The two observed departure rules generate a formal graph of 630 eight-node components.',
        limitations=['No new ODE integration','No infinite-time connection theorem',
            'Formal graph components need not be disconnected parts of an actual global dynamical network',
            'Symmetry copies of one saddle are not distinct newly discovered attractors'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=record['status'],nodes=graph['node_count'],edges=graph['edge_count'],
        components=graph['strongly_connected_component_count'],wall_seconds=record['wall_seconds'])))


if __name__ == '__main__':
    main()

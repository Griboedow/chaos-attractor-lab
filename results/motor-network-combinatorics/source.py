"""Exact permutation graph of the numerical motor saddle connections.

Group enumeration uses integer tuples only. The separate phase audit uses
the frozen numerical orbit and cannot certify isotropy or connections.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
IDENTITY = tuple(range(8))
PLUS = (5, 4, 6, 7, 2, 3, 0, 1)
MINUS = (5, 4, 7, 6, 2, 3, 1, 0)
HALF = (1, 0, 3, 2, 5, 4, 7, 6)


def compose(a, b):
    """Function composition a after b: (a*b)[j] = a[b[j]]."""
    return tuple(a[b[j]] for j in range(8))


def inverse(a):
    return tuple(a.index(j) for j in range(8))


def swap(a, b):
    p = list(IDENTITY)
    p[a], p[b] = p[b], p[a]
    return tuple(p)


def power(a, n):
    result = IDENTITY
    for _ in range(n):
        result = compose(a, result)
    return result


def closure(generators):
    group, pending = {IDENTITY}, [IDENTITY]
    while pending:
        a = pending.pop()
        for b in generators:
            p = compose(a, b)
            if p not in group:
                group.add(p)
                pending.append(p)
    assert all(tuple(sorted(p)) == IDENTITY for p in group)
    assert all(inverse(p) in group for p in group)
    return group


def pair_action(a):
    assert all(a[2*j] // 2 == a[2*j+1] // 2 for j in range(4))
    return tuple(a[2*j] // 2 for j in range(4))


def conjugate(a, b):
    return compose(compose(a, b), inverse(a))


def physical(states):
    x = np.concatenate((states[:, 16:], -states[:, 16:].sum(axis=1, keepdims=True)), axis=1)
    return np.stack((states[:, :8], states[:, 8:16], x), axis=2)


def exact_graph():
    h0 = closure((swap(2, 3), swap(6, 7)))
    h = closure((*h0, HALF))
    g = closure((*h, PLUS, MINUS))
    kernel = {a for a in g if pair_action(a) == tuple(range(4))}
    flips = closure(tuple(swap(2*j, 2*j+1) for j in range(4)))
    rotation = tuple(2*pair_action(PLUS)[j//2] + j % 2 for j in range(8))
    normalizer = {a for a in g if {conjugate(a, b) for b in h0} == h0}
    assert len(g) == 64 and len(h) == 8 and len(h0) == 4
    assert kernel == flips and len(kernel) == 16
    assert normalizer == kernel
    assert rotation in g and power(rotation, 4) == IDENTITY
    assert power(rotation, 2) != IDENTITY
    assert all(conjugate(a, b) in kernel for a in g for b in kernel)
    assert {compose(b, power(rotation, k)) for b in kernel for k in range(4)} == g
    assert all(compose(a, b) == compose(b, a) for a in kernel for b in kernel)
    assert all(power(a, 2) == IDENTITY for a in kernel)
    assert HALF not in h0 and power(HALF, 2) in h0
    cosets = {tuple(sorted(compose(b, a) for b in h)) for a in g}
    levels = [pair_action(power(rotation, k)) for k in range(4)]
    ordered = []
    for level, action in enumerate(levels):
        members = sorted(c for c in cosets if pair_action(c[0]) == action)
        assert len(members) == 2
        for index, c in enumerate(members):
            ordered.append((f"L{level}{'ab'[index]}", level, c))
    node_of = {p: name for name, _, c in ordered for p in c}
    level_of = {name: level for name, level, _ in ordered}
    nodes, edges = [], []
    for name, level, c in ordered:
        representative = c[0]
        targets = {node_of[compose(departure, representative)] for departure in (PLUS, MINUS)}
        for a in c:
            assert {node_of[compose(departure, a)] for departure in (PLUS, MINUS)} == targets
        assert len(targets) == 2
        assert all(level_of[target] == (level + 1) % 4 for target in targets)
        nodes.append(dict(id=name, level=level, representative=representative,
            left_coset=c, pair_action=pair_action(representative),
            expanding_pair=[j for j in range(8) if representative[j] in (6, 7)],
            contracting_pair=[j for j in range(8) if representative[j] in (2, 3)],
            outgoing=sorted(targets),
            representative_plus_target=node_of[compose(PLUS, representative)],
            representative_minus_target=node_of[compose(MINUS, representative)]))
        edges += [dict(source=name, target=target) for target in sorted(targets)]
    adjacency = {node['id']: node['outgoing'] for node in nodes}
    assert all(sum(edge['target'] == node['id'] for edge in edges) == 2 for node in nodes)
    cycles = set()

    def visit(path):
        for destination in adjacency[path[-1]]:
            if destination == path[0]:
                cycle = tuple(path)
                cycles.add(min(cycle[k:] + cycle[:k] for k in range(len(cycle))))
            elif destination not in path:
                visit(path + [destination])

    for node in nodes:
        visit([node['id']])
    assert sum(len(c) == 4 for c in cycles) == 16
    assert sum(len(c) == 8 for c in cycles) == 8
    selected = {}
    for name, generator in (('plus', PLUS), ('minus', MINUS)):
        sequence = [node_of[power(generator, k)] for k in range(5)]
        assert len(set(sequence[:4])) == 4 and sequence[0] == sequence[4]
        selected[name] = dict(permutation_order=next(k for k in range(1, 65)
            if power(generator, k) == IDENTITY), node_sequence=sequence,
            fourth_power=power(generator, 4))
    return dict(status='exact_integer_checks_passed',
        action_convention='(q dot a)[j]=q[a[j]]; (q dot a) dot b = q dot compose(a,b)',
        node_convention='Left cosets H a. Transported departure: H a -> H compose(g,a).',
        composition='compose(a,b)[j]=a[b[j]]',
        generators=dict(plus=PLUS, minus=MINUS, half=HALF),
        group=sorted(g), pointwise_subgroup=sorted(h0), known_orbit_subgroup=sorted(h),
        pair_flip_kernel=sorted(kernel), pointwise_normalizer_in_group=sorted(normalizer),
        orders=dict(group=len(g), pointwise_subgroup=len(h0), known_orbit_subgroup=len(h),
            kernel=len(kernel), normalizer=len(normalizer), nodes=len(nodes), edges=len(edges)),
        semidirect_product=dict(kernel='C2^4', quotient='C4', order_four_lift=rotation,
            lift_pair_action=pair_action(rotation), decomposition='kernel times lift powers is all 64 elements'),
        conditional_isotropy_lemma=dict(
            assumptions=['A nonconstant periodic orbit of the permutation-equivariant ODE exists.',
                'Its six cluster states are pairwise distinct, with pairs (2,3) and (6,7).',
                'The permutation h maps this orbit to itself with a time shift.'],
            normalizer_quotient_order=len(normalizer)//len(h0),
            normalizer_quotient_exponent=2,
            conclusion='The full orbital stabilizer intersected with G equals H; exactly eight distinct G-image periodic orbits. This does not determine the full S8 orbital stabilizer.'),
        nodes=nodes, edges=edges, simple_cycles=sorted(cycles),
        simple_cycle_counts={'4': 16, '8': 8}, selected_generator_circuits=selected,
        limitation='Edges encode observed finite-time destination permutations and their symmetry images. No infinite-time connection, network attraction, or chaotic-attractor theorem is supplied.')


def phase_audit(graph, source):
    with np.load(source) as data:
        times = data['time'].copy()
        rows = physical(data['full_states'])
    period = float(times[-1])
    original_closure = float(np.max(abs(rows[-1] - rows[0])))
    rows[-1] = rows[0]
    spline = CubicSpline(times, rows, bc_type='periodic')
    phase_grid = times[:-1]
    grid_rows = rows[:-1]
    h = {tuple(a) for a in graph['known_orbit_subgroup']}
    records = []
    for a in graph['group']:
        target = rows[0, list(a)]
        scores = np.sum((grid_rows - target)**2, axis=(1, 2))
        minima = np.flatnonzero((scores <= np.roll(scores, 1)) & (scores <= np.roll(scores, -1)))
        candidates = [(float(scores[k]), float(phase_grid[k])) for k in minima]

        def objective(phase):
            return float(np.sum((spline(phase % period) - target)**2))

        for k in minima:
            left = phase_grid[k-1] if k else phase_grid[-1] - period
            right = phase_grid[k+1] if k+1 < len(phase_grid) else period
            result = minimize_scalar(objective, bounds=(left, right), method='bounded',
                                     options={'xatol': 1e-13})
            candidates.append((float(result.fun), float(result.x % period)))
        score, phase = min(candidates)
        residual = spline(phase % period) - target
        records.append(dict(permutation=a, in_known_H=tuple(a) in h, phase=phase,
            normalized_phase=phase/period, euclidean_residual=float(np.sqrt(score)),
            max_coordinate_residual=float(np.max(abs(residual))),
            local_minima_refined=len(minima)))
    inside = [v['euclidean_residual'] for v in records if v['in_known_H']]
    outside = [v['euclidean_residual'] for v in records if not v['in_known_H']]
    assert len(inside) == 8 and len(outside) == 56
    assert max(inside) < 1e-8 and min(outside) > 1e-3
    return dict(status='numerical_phase_audit_passed', period=period, records=records,
        input_grid_points=len(times), original_endpoint_closure_max=original_closure,
        maximum_H_residual=max(inside), minimum_non_H_residual=min(outside),
        endpoint_forced_periodic_for_interpolation=True,
        method='For each of 64 permutations, refine every local minimum of the 1000-phase sampled squared distance using periodic cubic interpolation. Retain original grid candidates too.',
        limitation='Numerical interpolation and scalar minimization only; not an interval exclusion over all phases. No fresh ODE integration.')


def draw(graph, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    colors = ['#234D78', '#277D75', '#936431', '#76549A']
    centers = [(0, 0), (3.2, 0), (3.2, -2.5), (0, -2.5)]
    coords = {}
    for node in graph['nodes']:
        x, y = centers[node['level']]
        shift = -.38 if node['id'].endswith('a') else .38
        coords[node['id']] = (x + shift, y)
    fig, ax = plt.subplots(figsize=(9, 6))
    for edge in graph['edges']:
        start, end = coords[edge['source']], coords[edge['target']]
        arrow = FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=12,
            shrinkA=23, shrinkB=23, color='#78838F', alpha=.7, linewidth=1.2,
            connectionstyle='arc3,rad=0.04')
        ax.add_patch(arrow)
    for node in graph['nodes']:
        x, y = coords[node['id']]
        ax.scatter([x], [y], s=1800, color=colors[node['level']], zorder=3,
                   edgecolor='white', linewidth=1.7)
        ax.text(x, y, node['id'], color='white', ha='center', va='center', fontsize=11, zorder=4)
    for level, (x, y) in enumerate(centers):
        node = next(n for n in graph['nodes'] if n['level'] == level)
        label = 'Expanding pair ' + ','.join(str(v) for v in node['expanding_pair'])
        label += '\nContracting pair ' + ','.join(str(v) for v in node['contracting_pair'])
        ax.text(x, y + .58, label, ha='center', va='center', fontsize=10, color=colors[level])
    ax.set(xlim=(-1.25, 4.45), ylim=(-3.25, 1.1))
    ax.axis('off')
    fig.suptitle('Candidate symmetry graph: eight periodic nodes, two departures each', fontsize=13)
    fig.text(.5, .025, 'Exact permutation combinatorics; edge existence and network attraction are not proved.\nRow indices are zero-based. The a/b labels distinguish two orbit cosets at each level.',
             ha='center', fontsize=9, color='#414B56')
    fig.tight_layout(rect=(0, .08, 1, .95))
    fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/motor-network-combinatorics')
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError('Preserve evidence: choose a new output directory')
    started = time.monotonic()
    source = ROOT/'results/motor-postcritical/six-cluster-cycle-rho-3p268/cycle.npz'
    files = [source, ROOT/'results/motor-postcritical/saddle-escape/network-algebra.json',
             ROOT/'results/motor-postcritical/saddle-escape/summary.json']
    source_bytes = Path(__file__).read_bytes()
    graph = exact_graph()
    audit = phase_audit(graph, source)
    output.mkdir(parents=True)
    (output/'source.py').write_bytes(source_bytes)
    (output/'exact-graph.json').write_text(json.dumps(graph, indent=2)+'\n', encoding='utf-8')
    (output/'phase-audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    draw(graph, output/'candidate-graph.png')
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status='passed', wall_seconds=time.monotonic()-started, python=platform.python_version(),
        source_sha256=hashlib.sha256(source_bytes).hexdigest(),
        inputs={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
        outputs={path.name:hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in output.iterdir() if path.is_file()},
        group_order=64, known_H_order=8, graph_nodes=8, graph_edges=16,
        proof_scope='Exact finite permutation group and formal graph; conditional isotropy lemma. Phase audit is numerical. No new ODE simulation and no connection theorem.')
    (output/'analysis.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=record['status'], group_order=64, nodes=8, edges=16,
        phase_H_max=audit['maximum_H_residual'], phase_non_H_min=audit['minimum_non_H_residual'],
        wall_seconds=record['wall_seconds'])))


if __name__ == '__main__':
    main()

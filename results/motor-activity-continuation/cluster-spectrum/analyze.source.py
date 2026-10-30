"""Audit the completed 14+3+3+3 conditional spectrum and draw its evidence."""
from datetime import datetime, timezone
from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_determinant(matrix):
    a = [[Fraction(float(value)) for value in row] for row in matrix]
    determinant = Fraction(1)
    for j in range(len(a)):
        pivot = next(i for i in range(j, len(a)) if a[i][j])
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            determinant = -determinant
        value = a[j][j]
        determinant *= value
        for i in range(j + 1, len(a)):
            factor = a[i][j] / value
            for k in range(j, len(a)):
                a[i][k] -= factor * a[j][k]
    return determinant


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    folder = args.root.resolve()
    run = folder / 'lane-46'
    out = folder / 'analysis'
    if out.exists():
        raise FileExistsError('Preserve existing analysis')
    record = json.loads((run / 'analysis.json').read_text())
    parent = json.loads((run / 'parent-2000-analysis.json').read_text())
    checks = dict(complete=record['status'] == 'complete', requested_duration=record['observed_tau'] == 4000.,
                  warmup_duration=record['elapsed_tau'] - record['observed_tau'] == 200.,
                  eight_complete_windows=len(record['windows']) == 8 and all(w['complete'] and w['duration'] == 500. for w in record['windows']),
                  exact_groups=record['groups'] == [[0], [1, 2], [3, 6], [4, 7], [5]],
                  exact_weights=record['weights'] == [.125, .25, .25, .25, .125],
                  same_parent_projection=record['initial_projection_error'] == parent['initial_projection_error'],
                  parent_source_sha256=record['source_sha256'] == parent['source_sha256'],
                  parent_metadata_hash=sha(run / 'parent-2000-analysis.json') == record['parent_analysis_sha256'],
                  seed_hash=sha(run / 'seed.npz') == record['seed_sha256'],
                  history_hash=sha(run / 'history.npz') == record['history_sha256'],
                  decomposition=max(v for item in record['complete_decomposition_controls'] for v in item.values()) < 1e-10,
                  endpoint_decomposition=max(record['final_decomposition_control'].values()) < 1e-10,
                  short_full_map=record['independent_full_23D_short_map']['passed'],
                  trace_sum=max(record['block_trace_errors_per_tau']) < 1e-9,
                  local_trace=max(record['max_local_trace_errors']) < 1e-8,
                  local_flow=record['max_local_flow_error'] < 1e-8,
                  transported_flow=abs(record['flow_exponent'] - record['telescopic_flow_exponent']) < 1e-9,
                  weighted_volume=abs(record['weighted_volume_residual']) < 1e-6)
    for name, digest in record['source_hashes'].items():
        checks['source_' + name] = sha(run / name) == digest
    with np.load(run / 'seed.npz') as seed:
        determinant = exact_determinant(seed['full_tangent_transform'])
        checks['exact_basis_invertible'] = determinant == -512
    with np.load(run / 'history.npz') as data:
        times, states = data['time'], data['states']
        final_basis = data['final_quotient_basis']
        contrasts = data['final_contrast_bases']
        logs, traces = data['log_sums'], data['integrated_block_traces']
        checkpoints = data['checkpoints']
        checks['sampling_grid'] = bool(len(times) == 8400 and np.array_equal(times, np.arange(1, 8401) / 2))
        checks['all_finite'] = bool(np.isfinite(states).all() and np.isfinite(logs).all())
        checks['quotient_basis_orthogonal'] = bool(np.max(abs(final_basis.T @ final_basis - np.eye(14))) < 1e-12)
        checks['contrast_bases_orthogonal'] = bool(max(np.max(abs(q.T @ q - np.eye(3))) for q in contrasts) < 1e-12)
        checks['final_rates_match'] = bool(np.array_equal(logs / 4000, np.r_[record['quotient_spectrum'], np.asarray(record['contrast_spectra']).ravel()]))
        checks['determinant_trace'] = abs(float(logs.sum() - traces.sum())) < 1e-6
        checks['checkpoints_saved'] = checkpoints.shape == (9, 267)
        checkpoint_times = np.r_[200., 200. + 500. * np.arange(1, 9)]
        checks['checkpoint_times'] = bool(np.array_equal(checkpoints[:, 0], checkpoint_times))
        checks['checkpoint_states'] = bool(np.array_equal(checkpoints[:, 1:15], states[(2 * checkpoint_times).astype(int) - 1]))
        basis_errors = []
        for checkpoint in checkpoints:
            q = checkpoint[15:211].reshape(14, 14)
            basis_errors.append(np.max(abs(q.T @ q - np.eye(14))))
            for q in checkpoint[211:238].reshape(3, 3, 3):
                basis_errors.append(np.max(abs(q.T @ q - np.eye(3))))
        checks['all_checkpoint_bases_orthogonal'] = bool(max(basis_errors) < 1e-12)
    sample = states[times > 200]
    weights = np.array(record['weights'])
    positions = np.c_[sample[:, 10:], -(sample[:, 10:] @ weights[:-1]) / weights[-1]]
    pair_distances = []
    for i in range(5):
        for j in range(i + 1, 5):
            distances = np.sqrt((sample[:, i] - sample[:, j]) ** 2 + (sample[:, 5 + i] - sample[:, 5 + j]) ** 2 + (positions[:, i] - positions[:, j]) ** 2)
            pair_distances.append(dict(cluster_indices=[i, j], minimum=float(distances.min()),
                                       rms=float(np.sqrt(np.mean(distances ** 2))), maximum=float(distances.max())))
    checks['weighted_position_constraint'] = bool(np.max(abs(positions @ weights)) < 1e-12)
    assert all(checks.values()), {k: v for k, v in checks.items() if not v}
    quotient = np.array([w['quotient_spectrum'] for w in record['windows']])
    trans = np.array([w['contrast_spectra'] for w in record['windows']])
    t = np.arange(1, len(quotient) + 1) * 500
    spectrum = np.r_[record['quotient_spectrum'], np.asarray(record['contrast_spectra']).ravel()]
    out.mkdir(parents=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    axes[0].plot(t, quotient[:, 0], 'o-', label='Quotient QR 1', color='#A85636')
    axes[0].plot(t, quotient[:, 1], 'o-', label='Quotient QR 2', color='#747B85')
    for k, pair in enumerate(([1, 2], [3, 6], [4, 7])):
        axes[0].plot(t, trans[:, k, 0], '.-', label=f'Pair {pair[0]},{pair[1]} leading')
    axes[0].axhline(0, color='black', linewidth=.8)
    axes[0].set(xlabel='Observed tau (500-tau windows)', ylabel='Finite-window exponent',
                title='Expansion within clusters; transverse contraction')
    axes[0].legend(fontsize=8, ncol=2)
    axes[0].grid(alpha=.18)
    colors = ['#A85636'] * 14 + ['#477D99'] * 9
    axes[1].bar(np.arange(1, 24), spectrum, color=colors)
    axes[1].axhline(0, color='black', linewidth=.8)
    axes[1].axvline(14.5, color='#777777', linewidth=.8, linestyle='--')
    axes[1].set(xlabel='14 quotient + 9 pair-contrast directions', ylabel='Cumulative finite-time exponent',
                title='All 23 directions after 4,000 observed tau')
    axes[1].set_xticks([1, 7, 14, 15, 18, 21, 23])
    axes[1].grid(axis='y', alpha=.18)
    fig.suptitle('Lane 46: declared 1+2+2+2+1 cluster restriction', fontsize=12)
    fig.tight_layout()
    fig.savefig(out / 'spectrum.png', dpi=170)
    plt.close(fig)
    summary = dict(status='passed', created_utc=datetime.now(timezone.utc).isoformat(), checks=checks,
                   source_sha256=sha(Path(__file__)), run_record_sha256=sha(run / 'analysis.json'),
                   exact_transform_determinant=str(determinant), quotient_spectrum=record['quotient_spectrum'],
                   contrast_spectra=record['contrast_spectra'], sorted_union_spectrum=record['sorted_union_spectrum'],
                   minimum_window_quotient_first=float(quotient[:, 0].min()),
                   maximum_window_contrast_leaders=trans[:, :, 0].max(axis=0).tolist(),
                   different_cluster_pair_distances=pair_distances,
                   minimum_different_cluster_rms=min(item['rms'] for item in pair_distances),
                   flow_exponent=record['flow_exponent'], telescopic_flow_exponent=record['telescopic_flow_exponent'],
                   weighted_volume_residual=record['weighted_volume_residual'],
                   trace_errors_per_tau=record['block_trace_errors_per_tau'],
                   wall_seconds=record['wall_seconds'], figure_sha256=sha(out / 'spectrum.png'),
                   limitation='Finite numerical evidence for this declared restriction. No theorem of asymptotic chaos, transverse stability or full-space attraction.')
    write = lambda p, r: p.write_text(json.dumps(r, indent=2) + '\n', encoding='utf-8', newline='\n')
    write(out / 'analysis.json', summary)
    snapshot = folder / 'analyze.source.py'
    if snapshot.resolve() != Path(__file__).resolve():
        if snapshot.exists():
            raise FileExistsError('Do not overwrite an analysis source snapshot')
        shutil.copyfile(Path(__file__), snapshot)
    inventory = {p.relative_to(folder).as_posix(): dict(bytes=p.stat().st_size, sha256=sha(p))
                 for p in sorted(folder.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
    write(folder / 'manifest.json', dict(status='passed', created_utc=datetime.now(timezone.utc).isoformat(),
                                        files=inventory, file_count=len(inventory), checks=checks))
    print(json.dumps({k: summary[k] for k in ('status', 'quotient_spectrum', 'contrast_spectra', 'maximum_window_contrast_leaders')}))


if __name__ == '__main__':
    main()

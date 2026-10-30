"""Refine plotting distances without altering any integrated trajectory."""
import datetime
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import linear_sum_assignment, minimize_scalar

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'research'))
import motor_log_contrasts as log
import motor_pair_coordinates as pair


def history(folder):
    meta = json.loads((folder/'summary.json').read_text())
    if 'parent' not in meta:
        times, chart, norms, matches = pair.load_chain(folder)
        return times, chart, np.log10(norms), matches
    data = np.load(folder/'trajectory.npz')
    times, chart, logs, matches = history(ROOT/meta['parent'])
    own_matches = json.loads((folder/'diagnostics.json').read_text())['orbit_matches']
    return (np.r_[times, data['time'][1:]], np.r_[chart, log.to_chart(data['log_states'][1:])],
            np.r_[logs, data['log10_halfcontrast_norms'][1:]], matches+own_matches[1:])


def main():
    folder = HERE/'step-0p1-to-280k'
    output = HERE/'seventh-review'
    if (output/'audit.json').exists():
        raise FileExistsError(output/'audit.json')
    data = np.load(folder/'trajectory.npz')
    meta = json.loads((folder/'summary.json').read_text())
    diagnostic = json.loads((folder/'diagnostics.json').read_text())
    pt, ps, parent_logs, pm = history(ROOT/meta['parent'])
    times = np.r_[pt, data['time'][1:]]
    chart = np.r_[ps, log.to_chart(data['log_states'][1:])]
    logs = np.r_[parent_logs, data['log10_halfcontrast_norms'][1:]]
    old_matches = pm+diagnostic['orbit_matches'][1:]
    reference_path = ROOT/'results/motor-postcritical/six-cluster-cycle-rho-3p268/cycle.npz'
    reference = np.load(reference_path)
    ref_rows = pair.physical(reference['full_states'])
    ref_rows[-1] = ref_rows[0]
    period = float(reference['time'][-1])
    spline = CubicSpline(reference['time'], ref_rows, bc_type='periodic')
    grid = np.arange(1000)*period/1000
    grid_rows = spline(grid)
    dt = period/1000
    graph_path = ROOT/'results/motor-network-combinatorics/exact-graph.json'
    graph = json.loads(graph_path.read_text())
    lookup = {tuple(p): node['id'] for node in graph['nodes'] for p in node['left_coset']}
    matches, removed = [], []
    began = time.perf_counter()
    for old in old_matches:
        index = int(np.searchsorted(times, old['tau']))
        assert times[index] == old['tau']
        current = pair.chart_rows(chart[index])
        def cost(phase, details=False):
            rows = spline(phase % period)
            squared = np.sum((current[:, None, :]-rows[None, :, :])**2, axis=-1)
            ri, ci = linear_sum_assignment(squared)
            value = float(squared[ri, ci].sum())
            return (value, ri, ci) if details else value
        def newton(phase):
            for _ in range(8):
                value, ri, ci = cost(phase, True)
                rows, first, second = spline(phase % period), spline(phase % period, 1), spline(phase % period, 2)
                delta = current[ri]-rows[ci]
                denominator = float(np.sum(first[ci]**2)-np.sum(delta*second[ci]))
                if denominator <= 0:
                    break
                step = float(np.sum(delta*first[ci]))/denominator
                step = max(-dt, min(dt, step))
                if abs(step) < 1e-14:
                    break
                for _ in range(12):
                    if cost(phase+step) <= value:
                        phase += step
                        break
                    step *= .5
                else:
                    break
            return phase
        squares = np.sum((current[None, :, None, :]-grid_rows[:, None, :, :])**2, axis=-1)
        values = np.array([sq[linear_sum_assignment(sq)].sum() for sq in squares])
        minima = np.flatnonzero((values <= np.roll(values, 1)) & (values <= np.roll(values, -1)))
        guesses = [float(grid[np.argmin(values)]), float(old['phase'])]
        if matches:
            guesses.append((matches[-1]['phase']+old['tau']-matches[-1]['tau']) % period)
        for k in minima:
            fit = minimize_scalar(cost, bounds=(grid[k]-dt, grid[k]+dt), method='bounded', options={'xatol': 1e-13})
            guesses.append(float(fit.x))
        refined = [newton(phase) for phase in guesses]
        phase = min(refined, key=cost)
        squared, ri, ci = cost(phase, True)
        distance = math.sqrt(squared)
        assert distance <= old['distance']+1e-8
        node = lookup.get(tuple(ci))
        matches.append(dict(tau=old['tau'], distance=distance, phase=float(phase % period),
                            permutation_current_from_reference=ci.tolist(), candidate_graph_node=node,
                            original_distance=old['distance'], phase_grid_local_minima=len(minima)))
        if old['distance'] > 1e-3 and distance < 1e-6:
            removed.append(dict(tau=old['tau'], previous_distance=old['distance'], refined_distance=distance))
    visits = []
    for match in matches:
        if match['distance'] >= 1e-6 or match['candidate_graph_node'] is None:
            continue
        if not visits or visits[-1]['node'] != match['candidate_graph_node']:
            visits.append(dict(node=match['candidate_graph_node'], first_near_tau=match['tau'], last_near_tau=match['tau']))
        else:
            visits[-1]['last_near_tau'] = match['tau']
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='completed',
                  wall_seconds=time.perf_counter()-began, source_sha256=pair.sha(__file__),
                  trajectory_sha256=pair.sha(folder/'trajectory.npz'), original_diagnostic_sha256=pair.sha(folder/'diagnostics.json'),
                  reference_sha256=pair.sha(reference_path), graph_sha256=pair.sha(graph_path),
                  method='1000-phase grid; refine every sampled local minimum; include prior and propagated phase; safeguarded analytic Newton refinement with optimal row assignment.',
                  numerical_trajectory_unchanged=True, reference_interpolation_unchanged=True,
                  removed_large_matching_spikes=removed, visits=visits, matches=matches,
                  limitation='A numerical distance diagnostic, not a certified global phase minimum or heteroclinic proof.')
    pair.write(output/'audit.json', record)
    (output/'source.py').write_bytes(Path(__file__).read_bytes())
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    prediction = json.loads((HERE/'pre-event-prediction.json').read_text())
    departures = prediction['frozen_departure_times']+[r['recovery_upward_crossing_1e_minus_4'] for r in log.log_recoveries(folder)]
    gaps = np.diff(departures)
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), constrained_layout=True)
    for j, color in enumerate(['#4477AA', '#EE6677', '#228833', '#AA3377']):
        axes[0].plot(times, logs[:, j], color=color, lw=1.05, label=f'Rows {2*j},{2*j+1}')
    axes[0].axhline(math.log10(np.nextafter(0., 1.)), color='#777777', ls=':', label='Smallest positive FP64 value')
    axes[0].axvline(40000., color='#999999', ls='--', lw=.8)
    axes[0].axvline(150000., color='#999999', ls='--', lw=.8)
    axes[0].set(title='A contrast below the FP64 amplitude range returns and separates its pair', ylabel='log10 half-contrast norm', xlabel='Normalized time tau')
    axes[0].legend(ncol=3, loc='lower left')
    axes[1].semilogy([m['tau'] for m in matches], [m['distance'] for m in matches], color='#222222', lw=.9)
    axes[1].axhline(1e-9, color='#888888', ls=':', label='Approximate reference interpolation floor')
    axes[1].set(title=f'Refined phase matching: {len(departures)} finite-time departures', ylabel='Distance to permuted saddle', xlabel='Normalized time tau')
    axes[1].legend(loc='upper right')
    axes[2].plot(np.arange(2, len(gaps)+1), gaps[1:]/gaps[:-1], 'o-', color='#4477AA', label='Successive departure-gap ratio')
    axes[2].axhline(prediction['kappa'], color='#AA3377', ls='--', label='Saddle rate ratio |c| / e')
    axes[2].set(title='A conditional residence-growth diagnostic', xlabel='Later inter-departure gap index', ylabel='Ratio of consecutive gaps')
    axes[2].set_xticks(np.arange(2, len(gaps)+1))
    axes[2].legend()
    for ax in axes:
        ax.grid(alpha=.2)
    image = output/'log-contrast-diagnostic.png'
    fig.savefig(image, dpi=160)
    plt.close(fig)
    pair.write(output/'plot.json', dict(image_sha256=pair.sha(image), audit_sha256=pair.sha(output/'audit.json'),
                                       source_sha256=pair.sha(__file__), departure_times=departures,
                                       gaps=gaps.tolist(), gap_ratios=(gaps[1:]/gaps[:-1]).tolist()))
    print(json.dumps(dict(status=record['status'], wall_seconds=record['wall_seconds'],
                          removed_spikes=len(removed), visits=visits), indent=2))


if __name__ == '__main__':
    main()

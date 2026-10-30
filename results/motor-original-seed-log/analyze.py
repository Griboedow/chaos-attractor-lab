"""Compare two archived-start logarithmic runs without inferring a basin or LLE."""
import datetime
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys

import numpy as np

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
sys.path.insert(0, str(ROOT/'research'))
import motor_pair_coordinates as pair
import motor_log_contrasts as log
import motor_log_original_seed as experiment


def main():
    data, records = [], []
    for name in ['step-0p1', 'step-0p05']:
        folder = BASE/name
        d = np.load(folder/'trajectory.npz')
        s = json.loads((folder/'summary.json').read_text())
        assert s['status'] == 'complete' and s['trajectory_sha256'] == pair.sha(folder/'trajectory.npz')
        initial = np.load(BASE/'initial/trajectory.npz')
        assert d['initial_log_state'].tobytes() == initial['final_log_state'].tobytes()
        assert d['original_initial_state'].tobytes() == initial['original_initial_state'].tobytes()
        assert np.array_equal(d['parameters'], initial['parameters'])
        times, z, norms = d['time'], d['log_states'], d['log10_halfcontrast_norms']
        final_window = times >= times[-1]-1000
        reconstructed = pair.chart_rows(log.to_chart(z))[:, experiment.INVERSE]
        minimum_outside, equality_records = [], []
        protected = {tuple(sorted(x)) for x in experiment.ORIGINAL_PAIRS}
        for i in range(8):
            for j in range(i+1, 8):
                distances = np.linalg.norm(reconstructed[:, i]-reconstructed[:, j], axis=1)
                index = np.argmin(distances)
                equal = np.all(reconstructed[:, i] == reconstructed[:, j], axis=1)
                hits = np.flatnonzero(equal)
                if len(hits):
                    equality_records.append(dict(rows=[i,j], protected_by_log_chart=(i,j) in protected,
                                                 first_equal_checkpoint=float(times[hits[0]]),
                                                 equal_at_final_checkpoint=bool(equal[-1]),
                                                 equal_checkpoint_count=int(np.sum(equal))))
                if (i,j) not in protected:
                    minimum_outside.append(dict(rows=[i,j], distance=float(distances[index]), tau=float(times[index]),
                                               final_distance=float(distances[-1]),
                                               bitwise_equal_at_any_checkpoint=bool(len(hits)),
                                               bitwise_equal_at_final_checkpoint=bool(equal[-1])))
        per_pair = []
        for j in range(4):
            index = int(np.argmin(norms[:,j]))
            after = norms[index:,j]
            recovered = np.flatnonzero(after > -4.) if norms[index,j] < -10 else np.array([],dtype=int)
            episodes, armed, low, low_time = [], False, None, None
            for t, value in zip(times, norms[:,j]):
                if value < -10.:
                    armed = True
                    if low is None or value < low:
                        low, low_time = float(value), float(t)
                if armed and value > -4.:
                    episodes.append(dict(minimum_log10_norm=low, minimum_tau=low_time,
                                          recovery_checkpoint_above1e_minus4=float(t)))
                    armed, low, low_time = False, None, None
            per_pair.append(dict(original_pair=experiment.ORIGINAL_PAIRS[j], minimum_log10_norm=float(norms[index,j]),
                                  minimum_tau=float(times[index]), final_log10_norm=float(norms[-1,j]),
                                  final1000_log10_range=[float(np.min(norms[final_window,j])),float(np.max(norms[final_window,j]))],
                                  first_postminimum_checkpoint_above1e_minus4=float(times[index+recovered[0]]) if len(recovered) else None,
                                  sampled_deep_contraction_recovery_episodes=episodes))
        bounds=[]
        for j,item in enumerate(s['underflow']):
            if item['maximum_omitted_component_log_bound'] is not None:
                with localcontext() as context:
                    context.prec=60
                    bound=Decimal(str(item['maximum_omitted_component_log_bound'])).exp()
                bounds.append(dict(original_pair=experiment.ORIGINAL_PAIRS[j],
                                   maximum_recorded_log_bound=item['maximum_omitted_component_log_bound'],
                                   corresponding_decimal=str(bound)))
        records.append(dict(run=name, parameters=s['parameters'], status=s['status'], final_tau=float(times[-1]),
                            wall_seconds=s['wall_seconds'], initial27_bytes_unchanged=True,
                            maximum_gauge_norm_squared_drift=s['maximum_gauge_norm_squared_drift'],
                            underflow=s['underflow'], omitted_mean_component_bounds=bounds, pairs=per_pair,
                            final1000_maximum_mean_coordinate_range=float(np.max(np.ptp(z[final_window,:11],axis=0))),
                            final1000_maximum_log_norm_range=float(np.max(np.ptp(norms[final_window],axis=0))),
                            all_four_below1e_minus4_at_final=bool(np.all(norms[-1] < -4.)),
                            final1000_all_four_below1e_minus4_checkpoint_count=int(np.sum(np.all(norms[final_window]<-4.,axis=1))),
                            minimum_sampled_unprotected_pair_separation=min(minimum_outside,key=lambda item:item['distance']),
                            all_unprotected_pair_separations=minimum_outside,
                            reconstructed_bitwise_equal_pairs=equality_records,
                            equality_scope='Reconstructed Cartesian rows only. Equality of protected pairs can coexist with nonzero logarithmic contrasts. Equality of unprotected pairs identifies precision limitations outside this chart.'))
        data.append(d)
    a,b=data
    assert np.array_equal(a['time'],b['time'])
    comparisons=[]
    for end in [20.,50.,100.,200.,500.,1000.,5000.,10000.]:
        if end > a['time'][-1]:
            continue
        mask=a['time']<=end
        comparisons.append(dict(through_tau=end,
                                 maximum_mean_difference=float(np.max(abs(a['log_states'][mask,:11]-b['log_states'][mask,:11]))),
                                 maximum_log10_norm_difference=float(np.max(abs(a['log10_halfcontrast_norms'][mask]-b['log10_halfcontrast_norms'][mask])))))
    result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='completed',
                runs=records, pointwise_controls=comparisons, source_sha256=pair.sha(__file__),
                initial_audit_sha256=pair.sha(BASE/'initial/initial-audit.json'),
                limitations=['Fixed four selected pair contrasts only.', 'Floating conversion error is documented; distant trajectories may separate.',
                             'No basin estimate, Lyapunov estimate, attracting-set theorem, or restoration of previously lost digits.'])
    pair.write(BASE/'comparison.json',result)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(11,7),sharex=True,constrained_layout=True)
    colors=['#4477AA','#EE6677','#228833','#AA3377']
    for ax,d,record in zip(axes,data,records):
        for j,color in enumerate(colors):
            ax.plot(d['time'],d['log10_halfcontrast_norms'][:,j],lw=1,color=color,label=f'Original rows {experiment.ORIGINAL_PAIRS[j]}')
        ax.set_title(record['run']+'; same archived row-46 initial state')
        ax.set_ylabel('log10 half-contrast norm'); ax.grid(alpha=.2); ax.legend(ncol=2,fontsize=9)
    axes[-1].set_xlabel('Normalized time tau')
    fig.savefig(BASE/'original-seed-log.png',dpi=160)
    plt.close(fig)
    entries=[]
    for path in sorted(BASE.rglob('*')):
        if path.is_file() and path.name!='manifest.json' and '__pycache__' not in path.parts:
            entries.append(dict(path=path.relative_to(BASE).as_posix(),sha256=pair.sha(path),bytes=path.stat().st_size))
    pair.write(BASE/'manifest.json',dict(created_utc=result['created_utc'],files=entries,file_count=len(entries),total_bytes=sum(x['bytes'] for x in entries)))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()

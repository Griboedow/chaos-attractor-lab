"""Audit and summarize the immutable CPU continuation once it has stopped."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def pair_statistics(states):
    x=np.c_[states[:,16:],-states[:,16:].sum(axis=1)]
    result=[]
    for i in range(8):
        for j in range(i+1,8):
            distances=np.sqrt((states[:,i]-states[:,j])**2+(states[:,8+i]-states[:,8+j])**2+(x[:,i]-x[:,j])**2)
            result.append(dict(rows=[i,j],minimum=float(distances.min()),maximum=float(distances.max()),
                rms=float(np.sqrt(np.mean(distances**2)))))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,default=Path(__file__).parent/'run-01')
    parser.add_argument('--output',type=Path,default=Path(__file__).parent/'analysis')
    args=parser.parse_args();run=args.run.resolve();out=args.output.resolve()
    if out.exists(): raise FileExistsError('Preserve previous analysis')
    report=json.loads((run/'analysis.json').read_text())
    old=json.loads((run/'input.json').read_text())
    assert report['status'] in ('complete','wall_or_absolute_cutoff','near_equilibrium')
    controls={}
    controls['source_script_snapshot']=all(sha(run/name)==value for name,value in report['source_hashes'].items())
    controls['input_metadata_snapshot']=sha(run/'input.json')==report['input_hashes']['results/motor-rows/activity-path/cpu-N8-row69-extended.json']
    controls['input_npz_snapshot']=sha(run/'input.npz')==report['input_hashes']['results/motor-rows/activity-path/cpu-N8-row69-extended.npz']
    controls['history_windows_unchanged']=report['windows'][:len(old['windows'])]==old['windows']
    controls['history_npz_prefixes']=all(report['inherited_history_prefix_bitwise_checks'].values())
    controls['final_endpoint_algebra']=report['final_endpoint_control']['passed']
    controls['final_npz_hash']=sha(run/'continuation.npz')==report['continuation_npz_sha256']
    controls['QR_trace_consistency']=report['max_local_trace_error']<1e-8 and report['sum_trace_error']<1e-9
    controls['transported_flow_consistency']=report['max_local_flow_identity_error']<1e-8
    controls['flow_log_telescoping']=abs(report['flow_exponent']-report['telescopic_flow_exponent'])<1e-9
    controls['restored_observation_flow_origin']=report['restoration']['recovered_flow_rate_error']<1e-12
    blocks=[]
    for i in range(4):
        windows=old['windows'][10*i:10*(i+1)]
        weight=np.array([v['duration'] for v in windows])
        assert abs(weight.sum()-1000)<1e-9
        rates=np.array([v['spectrum'] for v in windows])
        earlier=old['windows'][:10*(i+1)]
        blocks.append(dict(source='inherited',observed_end=(i+1)*1000,duration=1000.,
            spectrum=np.average(rates,axis=0,weights=weight).tolist(),
            cumulative_spectrum=np.average([v['spectrum'] for v in earlier],axis=0,
                weights=[v['duration'] for v in earlier]).tolist()))
    for item in report['new_blocks']:
        path=run/f"block-{item['index']:02d}.npz"
        controls[f"block_{item['index']}_npz_hash"]=sha(path)==item['npz_sha256']
        with np.load(path) as data:
            states=data['states'].copy()
        pairs=pair_statistics(states)
        blocks.append(dict(source='continued',observed_end=item['end_observed_tau'],duration=item['duration'],
            spectrum=item['spectrum'],cumulative_spectrum=item['cumulative_spectrum'],
            pair_distances=pairs,minimum_pair_rms=min(v['rms'] for v in pairs),
            smallest_instantaneous_pair_distance=min(v['minimum'] for v in pairs),
            flow_exponent=item['flow_exponent'],telescopic_flow_exponent=item['telescopic_flow_exponent'],
            weighted_volume_residual=item['weighted_volume_residual']))
    with np.load(run/'continuation.npz') as data:
        final=data['final_state'].copy();basis=data['final_basis'].copy()
        states=data['states'].copy();times=data['time'].copy()
        returns=data['section_states'].copy()
        inherited_arrays=dict(time=times,states=states,section_states=returns,
            section_time=data['section_time'].copy(),checkpoint=data['checkpoint'].copy())
    with np.load(run/'input.npz') as inherited:
        for name,values in inherited_arrays.items():
            original=inherited[name]
            controls['independent_bitwise_prefix_'+name]=values[:len(original)].tobytes()==original.tobytes()
    rows=np.column_stack((final[:8],final[8:16],np.r_[final[16:],-final[16:].sum()]))
    equal=[(i,j) for i in range(8) for j in range(i+1,8) if np.array_equal(rows[i],rows[j])]
    new=[v for v in blocks if v['source']=='continued']
    tail=returns[-min(1000,len(returns)):]
    return_rows=np.c_[tail,-tail[:,16:].sum(axis=1)]
    recurrence=[]
    for lag in range(1,min(65,len(return_rows)//2)):
        distances=np.linalg.norm(return_rows[lag:]-return_rows[:-lag],axis=1)
        recurrence.append(dict(lag=lag,minimum=float(distances.min()),median=float(np.median(distances)),
            maximum=float(distances.max()),rms=float(np.sqrt(np.mean(distances**2)))))
    spectrum=np.asarray(report['spectrum']);cumulative=np.array([b['cumulative_spectrum'] for b in blocks])
    controls['all_finite']=bool(np.isfinite(states).all() and np.isfinite(spectrum).all())
    controls['basis_orthogonality']=bool(np.max(abs(basis.T@basis-np.eye(23)))<1e-12)
    controls['weighted_volume_check']=abs(report['new_segment']['identity_residual'])<1e-6
    assert all(controls.values()),controls
    out.mkdir(parents=True)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),gridspec_kw={'width_ratios':[1.55,1]})
    tau=[b['observed_end'] for b in blocks]
    for j in range(8):
        axes[0].plot(tau,cumulative[:,j],label=f'QR {j+1}',linewidth=1.5)
    axes[0].axhline(0,color='black',lw=.8)
    axes[0].axvline(4000,color='#777777',lw=.8,ls='--')
    axes[0].set(xlabel='Total observed original tau',ylabel='Cumulative finite-time exponent',
        title='Exact endpoint/basis continuation from 4,000 tau')
    axes[0].legend(ncol=4,fontsize=8);axes[0].grid(alpha=.18)
    colors=['#AC5937' if value>.01 else '#6B7680' for value in spectrum]
    axes[1].bar(np.arange(1,24),spectrum,color=colors)
    axes[1].axhline(0,color='black',lw=.8)
    axes[1].set(xlabel='QR direction',ylabel='Finite-time exponent',
        title=f"All 23 directions after {report['observed_tau']:,.0f} tau")
    axes[1].set_xticks([1,6,7,12,18,23]);axes[1].grid(axis='y',alpha=.18)
    fig.suptitle('Published eight-row model: one retained trajectory, no tangent reset',fontsize=12)
    fig.tight_layout();fig.savefig(out/'cpu-spectrum.png',dpi=170);plt.close(fig)
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='passed',controls=controls,
        input_report_sha256=sha(run/'analysis.json'),input_npz_sha256=sha(run/'continuation.npz'),
        source_sha256=sha(Path(__file__)),source_status=report['status'],observed_tau=report['observed_tau'],
        added_tau=report['added_tau'],wall_seconds=report['wall_seconds'],spectrum=report['spectrum'],
        initial_4000_spectrum=old['spectrum'],thousand_tau_blocks=blocks,
        minimum_new_sixth_window_exponent=min(v['spectrum'][5] for v in new),
        maximum_new_seventh_window_exponent=max(v['spectrum'][6] for v in new),
        final_bitwise_equal_row_pairs=equal,
        final_basis_orthogonality_max=float(np.max(abs(basis.T@basis-np.eye(23)))),
        minimum_new_window_pair_rms=min(v['minimum_pair_rms'] for v in new),
        final_window_minimum_pair_rms=new[-1]['minimum_pair_rms'],
        poincare_recurrence=dict(tail_returns=len(tail),lags=recurrence,
            best_rms_lag=min(recurrence,key=lambda v:v['rms']),
            scope='Full labeled physical-row distance over the last up-to-1000 section returns, for lags 1 through 64. A finite numerical short-repetition diagnostic, not an aperiodicity theorem.'),
        flow_exponent=report['flow_exponent'],telescopic_flow_exponent=report['telescopic_flow_exponent'],
        trace_qr_error_per_tau=report['sum_trace_error'],new_segment=report['new_segment'],
        max_local_trace_error=report['max_local_trace_error'],
        max_local_flow_identity_error=report['max_local_flow_identity_error'],
        max_new_local_weighted_volume_error=report['max_new_weighted_volume_error'],
        sampled_domain=report['sampled_domain'],figure_sha256=sha(out/'cpu-spectrum.png'),
        limitation='Finite-time numerical continuation and consistency controls. A positive exponent estimate is not an asymptotic chaos or attraction theorem.')
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ('status','observed_tau','spectrum','minimum_new_sixth_window_exponent',
        'final_window_minimum_pair_rms','flow_exponent','trace_qr_error_per_tau')}))


if __name__=='__main__': main()

"""Test short return periods on a separate continuation of a cluster candidate."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import motor_postcritical_cycles as cycles


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output')
    meta=json.loads((args.input/'analysis.json').read_text())
    source=np.load(args.input/'history.npz')
    groups=meta['groups'];w=np.array(meta['weights']);p=np.array(meta['parameters'])
    initial=source['final_state'].copy();k=len(w)
    def section(t,y):
        return y[2*k]
    section.direction=1
    sol=solve_ivp(lambda t,y:cycles.field_jac(y,p,w)[0],(0.,1000.),initial,
        method='DOP853',rtol=2e-12,atol=2e-14,max_step=.025,events=section,dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    times=sol.t_events[0];points=sol.y_events[0]
    full=np.array([cycles.lift(y,groups) for y in points])
    lags=[]
    for lag in range(1,min(64,len(times)-1)+1):
        differences=np.linalg.norm(full[lag:]-full[:-lag],axis=1)
        lags.append(dict(lag=lag,comparisons=len(differences),minimum_full_reduced_distance=float(differences.min()),
            rms_full_reduced_distance=float(np.sqrt(np.mean(differences**2))),
            final_return_max_error=float(np.max(abs(full[-1]-full[-1-lag]))),
            final_period_guess=float(times[-1]-times[-1-lag])))
    intervals=np.diff(times)
    args.output.mkdir(parents=True)
    sample_times=np.linspace(0.,1000.,10001)
    np.savez_compressed(args.output/'sample.npz',parameters=p,weights=w,time=sample_times,
        cluster_states=sol.sol(sample_times).T,section_time=times,section_states=full,
        cluster_section_states=points,initial_state=initial,final_state=sol.y[:,-1])
    report=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='complete',
        input_hashes={name:sha(args.input/name) for name in ['analysis.json','history.npz']},
        source_hashes={Path(__file__).name:sha(Path(__file__)),Path(cycles.__file__).name:sha(Path(cycles.__file__))},
        input_folder=str(args.input),groups=groups,weights=w.tolist(),parameters=p.tolist(),duration=1000.,
        exact_saved_endpoint_used=bool(np.array_equal(initial,source['final_state'])),
        section='Positive crossing of the first cluster displacement through zero; all roots retained.',
        crossing_count=len(times),interval_min=float(intervals.min()),interval_max=float(intervals.max()),
        interval_std=float(intervals.std()),lags=lags,
        sample_sha256=sha(args.output/'sample.npz'),
        limitations=['Finite return search at lags1 through64; no proof of aperiodicity or exclusion of longer periods.',
            'A separate tighter state continuation, not an extension of the archived QR accumulation.',
            'No interpolation-based orbit closure is treated as a certificate.'])
    (args.output/'analysis.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    for src in [Path(__file__),Path(cycles.__file__)]:
        (args.output/src.name).write_bytes(src.read_bytes())
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    axes[0].plot(times[1:],intervals,'.',ms=3)
    axes[0].set(xlabel='Time from archived spectrum endpoint',ylabel='Consecutive section-return interval',title='Every positive displacement crossing')
    axes[1].scatter(intervals[:-1],intervals[1:],s=10,alpha=.65)
    axes[1].set(xlabel='Return interval n',ylabel='Return interval n+1',title='Finite return diagnostic; no attractor theorem')
    fig.suptitle('Five-cluster lane46: a separate1000-tau return check')
    fig.savefig(args.output/'return-diagnostic.png',dpi=160);plt.close(fig)
    print('complete',len(times),'crossings','best lag minimum',min(x['minimum_full_reduced_distance'] for x in lags))


if __name__=='__main__':
    main()

"""Audit every declared full-space probe without selecting successful outcomes."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

import numpy as np
from report_motor_activity_gpu import row_summary


ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'results/motor-activity-continuation/full-space-perturbations')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new audit folder')
    folder=args.input
    protocol=json.loads((folder/'protocol.json').read_text())
    metadata=json.loads((folder/'outcomes.json').read_text())
    data=np.load(folder/'outcomes.npz')
    baseline=ROOT/protocol['source_path']
    original=np.load(baseline)
    checks=[]
    def check(label,passed,**details):
        checks.append(dict(check=label,passed=bool(passed),**details))
    check('source_checkpoint_hash',sha(baseline)==protocol['source_sha256'])
    check('all20declared_cases',len(protocol['cases'])==20 and len(data['status'])==20)
    check('all20complete',np.all(data['status']==0))
    check('all20full_duration',np.all(abs(data['elapsed_original']-5000.)<1e-8))
    check('settings',metadata['settings']['time_mode']=='original' and metadata['settings']['transient_scaled']==1000. and metadata['settings']['observation_scaled']==4000.)
    for name,digest in protocol['source_snapshot_sha256'].items():
        check('executed_source_hash',sha(folder/name)==digest,source=name)
    outcomes=[]
    for i,case in enumerate(protocol['cases']):
        sid=case['source_lane'];j=int(np.flatnonzero(original['source_lane_indices']==sid)[0])
        start=original['final_states'][j]
        expected=start+np.asarray(case['exact_reduced_delta'])
        check('declared_initial_state',np.array_equal(data['initial_states'][i],expected),case=i)
        check('unchanged_parameters',np.array_equal(data['parameters'][i],original['parameters'][j]),case=i)
        delta=data['initial_states'][i]-start
        norm=float(np.linalg.norm(np.r_[delta,-sum(delta[16:])]))
        check('measured_physical_norm',norm==case['actual_full_physical_norm'],case=i)
        if case['sign']:
            check('requested_physical_norm',np.isclose(norm,case['requested_full_physical_norm'],rtol=1e-5,atol=0),case=i)
        else:
            check('baseline_state_bitwise',np.array_equal(data['initial_states'][i],start),case=i)
        outcomes.append(dict(case=i,source_lane=sid,sign=case['sign'],requested_norm=case['requested_full_physical_norm'],actual_norm=norm,
            mean_ftle=float(data['lle'][i]),four1000_tau_windows=data['lle_windows'][i].tolist(),
            initial_rows=row_summary(data['initial_states'][i]),final_rows=row_summary(data['final_states'][i])))
    args.output.mkdir(parents=True)
    report=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,outcomes=outcomes,
        source_sha256=sha(Path(__file__)),helper_sha256=sha(ROOT/'research/report_motor_activity_gpu.py'),
        input_sha256={p.name:sha(p) for p in [folder/'protocol.json',folder/'outcomes.json',folder/'outcomes.npz']},
        limitations=['All20declared cases retained; no random basin sampling.',
            'Finite-time positive growth can occur on periodic saddles and does not alone establish chaos.',
            'Returning near synchronized pairs does not certify asymptotic transverse stability.'])
    (args.output/'analysis.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    (args.output/'source.py').write_bytes(Path(__file__).read_bytes())
    (args.output/'row-summary-source.py').write_bytes((ROOT/'research/report_motor_activity_gpu.py').read_bytes())
    if report['status']!='passed':
        raise ArithmeticError([c for c in checks if not c['passed']])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    labels=['0','+1e-7','-1e-7','+1e-9','-1e-9']
    for ax,sid in zip(axes.flat,[36,42,46,49]):
        selected=[a for a in outcomes if a['source_lane']==sid]
        values=np.array([a['four1000_tau_windows'] for a in selected])
        for k in range(4):
            ax.scatter(np.arange(5)+.07*(k-1.5),values[:,k],s=22,label=f'Window{k+1}')
        ax.plot(np.arange(5),[a['mean_ftle'] for a in selected],color='black',marker='_',lw=1,label='Mean')
        ax.axhline(0,color='gray',lw=.7)
        ax.set(title=f'Original source lane {sid}',xticks=np.arange(5),xticklabels=labels,
            xlabel='Declared physical perturbation norm and sign',ylabel='Finite-time tangent growth per tau')
        ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8,ncol=3)
    fig.suptitle('Full-space perturbation probes: all20 outcomes\n1000 tau warmup +4000 observed; a positive rate alone does not identify an attractor')
    image=args.output/'perturbation-outcomes.png';fig.savefig(image,dpi=160);plt.close(fig)
    (args.output/'plot.json').write_text(json.dumps(dict(image_sha256=sha(image),analysis_sha256=sha(args.output/'analysis.json')),indent=2)+'\n',encoding='utf-8')
    print(report['status'],len(checks),'checks')


if __name__=='__main__':
    main()

"""Audit and plot the two completed continuations of the original fixed cohort."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import numpy as np


ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def row_summary(state):
    rows=np.stack([state[:8],state[8:16],np.r_[state[16:],-sum(state[16:])]],axis=1)
    pairs=[(i,j) for i in range(8) for j in range(i+1,8)]
    distances=np.array([np.linalg.norm(rows[i]-rows[j]) for i,j in pairs])
    exact=[list(p) for p in pairs if np.array_equal(rows[p[0]],rows[p[1]])]
    near=[list(p) for p,d in zip(pairs,distances) if d<1e-10]
    return dict(minimum_pair_distance=float(distances.min()),exact_equal_pairs=exact,near_pairs_below_1e_minus_10=near)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT/'results/motor-activity-continuation')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose a new output')
    baseline=ROOT/'results/motor-rows/activity-path/N8-lifetime-history.npz'
    old=np.load(baseline);ids=old['source_lane_indices'];checks=[];runs={}
    def check(label,passed,**details):checks.append(dict(check=label,passed=bool(passed),**details))
    for label in ['main','strict']:
        folder=args.root/('gpu-'+label);meta=json.loads((folder/'summary.json').read_text());data=np.load(folder/'history.npz')
        check('complete40blocks',meta['status']=='complete' and len(meta['blocks'])==40,setting=label)
        check('all32statescomplete',data['status'].shape==(40,32) and np.all(data['status']==0),setting=label)
        check('expectedwindowshape',data['lle_windows'].shape==(40,32,4),setting=label)
        check('samefixedcohort',np.array_equal(ids,data['source_lane_indices']),setting=label)
        check('sameparameters',np.array_equal(old['parameters'],data['parameters']),setting=label)
        check('inputhash',meta['source_sha256']==sha(baseline),setting=label)
        initial=np.load(folder/'initial.npz')
        check('exactinitialstates',np.array_equal(initial['states'],old['final_states']),setting=label)
        check('exactinitialtangents',np.array_equal(initial['tangents'],old['final_tangents']),setting=label)
        previous_states=old['final_states'];previous_tangents=old['final_tangents'];checkpoint_rows=[]
        for block in range(1,41):
            saved=np.load(folder/f'block-{block:02d}.npz')
            check('exactblockstatecarry',np.array_equal(saved['initial_states'],previous_states),setting=label,block=block)
            check('exactblocktangentinput',np.array_equal(saved['carried_initial_tangents'],previous_tangents),setting=label,block=block)
            check('completedblockduration',np.all(saved['elapsed_original']>=500.-1e-9),setting=label,block=block)
            previous_states=saved['final_states'];previous_tangents=saved['final_tangents']
            checkpoint_rows.append([row_summary(state) for state in previous_states])
        check('finalstatecheckpoint',np.array_equal(previous_states,data['final_states']),setting=label)
        check('finaltangentcheckpoint',np.array_equal(previous_tangents,data['final_tangents']),setting=label)
        for name,expected in meta['source_snapshot_sha256'].items():check('sourcearchive',sha(folder/name)==expected,setting=label,source=name)
        retained=np.all(data['lle']>.005,axis=0)&np.all(data['lle_windows']>.001,axis=(0,2))
        runs[label]=dict(folder=folder,data=data,meta=meta,retained=retained,checkpoint_rows=checkpoint_rows)
    lanes=[]
    for i,sid in enumerate(ids):
        lane=dict(source_lane=int(sid),parameters=old['parameters'][i].tolist(),previous5000_mean_ftle=float(old['lle'][:,i].mean()),initial_rows=row_summary(old['final_states'][i]))
        for label,r in runs.items():
            a=r['data'];lane[label]=dict(mean_over_new20000=float(a['lle'][:,i].mean()),minimum125window=float(a['lle_windows'][:,i].min()),
                maximum125window=float(a['lle_windows'][:,i].max()),retained_every_newwindow=bool(r['retained'][i]),
                final_rows=row_summary(a['final_states'][i]),
                saved_endpoint_rows=[dict(absolute_tau=7100.+500.*(k+1),**rows[i]) for k,rows in enumerate(r['checkpoint_rows'])])
        lanes.append(lane)
    report=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
        checks=checks,lanes=lanes,cohort_size=32,prior_absolute_tau=7100.,additional_tau=20000.,final_absolute_tau=27100.,
        block_tau=500.,window_tau=125.,windows_per_lane_per_setting=160,
        retained_counts={k:int(v['retained'].sum()) for k,v in runs.items()},
        retained_in_both_count=int(np.sum(runs['main']['retained']&runs['strict']['retained'])),
        source_sha256=sha(Path(__file__)),baseline_sha256=sha(baseline),
        input_sha256={str(p.relative_to(ROOT)):sha(p) for v in runs.values() for p in [v['folder']/'summary.json',v['folder']/'history.npz']},
        caveats=['A fixed selected cohort, not basin probabilities or independent random samples.',
            'Ten positive lanes already contain synchronized pairs; conditional transverse audits are separate.',
            'The stricter lane60 also develops synchronized pairs during this continuation; endpoint classes can change without losing positive tangent growth.',
            'The lone subthreshold lane49 is a numerically checked periodic saddle with three unstable directions, not a stable periodic capture.',
            'The settings compare finite-time persistence and distributions; chaotic endpoints need not agree pointwise.',
            'No asymptotic attractor or scientific-priority claim.'])
    args.output.mkdir(parents=True)
    (args.output/'analysis.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    (args.output/'source.py').write_bytes(Path(__file__).read_bytes())
    if report['status']!='passed':raise ArithmeticError([c for c in checks if not c['passed']])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(14,11),layout='constrained')
    fig.suptitle('Original eight-row activity cohort: two further 20,000-unit histories\nPositive tangent growth is a diagnostic; synchronized invariant sets require transverse checks',fontsize=13)
    for ax,label in zip(axes[0],['main','strict']):
        a=runs[label]['data'];im=ax.imshow(a['lle'].T,origin='lower',aspect='auto',cmap='viridis',vmin=-.01,vmax=.4,
            extent=[7100,27100,-.5,31.5]);fig.colorbar(im,ax=ax,label='Mean leading FTLE per 500 tau')
        ax.set(title=label.capitalize()+' settings; all 32 lanes retained in the archive',xlabel='Absolute normalized time tau',ylabel='Original source lane',
            yticks=np.arange(32),yticklabels=ids);ax.tick_params(axis='y',labelsize=7)
    ax=axes[1,0]
    for label,color,marker in [('main','#225ea8','o'),('strict','#d95f0e','x')]:
        ax.scatter(old['parameters'][:,0],runs[label]['data']['lle'].mean(axis=0),s=27,c=color,marker=marker,label=label)
    ax.axhline(0,color='black',lw=.7);ax.set(xlabel='rho (fixed source activity path)',ylabel='Leading FTLE over new 20,000 tau',title='Two settings; differences are not trajectory error bars');ax.legend();ax.grid(alpha=.2)
    ax=axes[1,1];x=np.arange(32)
    for offset,label,color in [(-.2,'main','#225ea8'),(.2,'strict','#d95f0e')]:
        counts=[len(l[label]['final_rows']['near_pairs_below_1e_minus_10']) for l in lanes]
        ax.bar(x+offset,counts,width=.36,label=label,color=color)
    ax.set(xticks=x,xticklabels=ids,xlabel='Original source lane',ylabel='Row pairs closer than 1e-10',title='Endpoint synchronization; no proof of transverse attraction')
    ax.tick_params(axis='x',labelsize=7,rotation=90);ax.legend();ax.grid(axis='y',alpha=.2)
    path=args.output/'cohort-continuation.png';fig.savefig(path,dpi=160);plt.close(fig)
    (args.output/'plot.json').write_text(json.dumps(dict(image_sha256=sha(path),analysis_sha256=sha(args.output/'analysis.json'),source_sha256=sha(Path(__file__))),indent=2)+'\n',encoding='utf-8')
    print(report['status'],len(checks),'checks',report['retained_counts'])


if __name__=='__main__':main()

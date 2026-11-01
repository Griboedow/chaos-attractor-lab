"""Freeze the completed three-row motor flow-quotient follow-up."""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    folders=[('full-spectrum','motor-three-two-positive-long'),('quotient-60000','motor-three-transverse-long'),('quotient-300000','motor-three-transverse-300k')]
    for label,name in folders:shutil.copytree(args.external/name,out/label)
    results=[]
    for label,_ in folders[1:]:
        source=out/label;protocol=json.loads((source/'protocol.json').read_text());a=json.loads((source/'analysis.json').read_text())
        original=Path(protocol['source']);assert sha(original)==protocol['sha256']
        h=np.load(source/'history.npz');windows=np.array([w['spectrum'] for w in a['windows']]);duration=protocol['duration']
        assert len(h['time'])==a['samples'] and np.allclose(h['log_sums']/duration,a['spectrum'],atol=1e-14,rtol=0)
        assert np.allclose(windows.mean(0),a['spectrum'],atol=1e-13,rtol=0)
        residual=abs(h['log_sums'].sum()+np.log(a['final_speed']/a['initial_speed'])-h['integrated_trace'])/duration
        assert abs(residual-a['divergence_identity_residual_per_tau'])<1e-15
        results.append(dict(run=label,duration=duration,leading=a['spectrum'][:2],second_window_min=float(windows[:,1].min()),
            second_window_max=float(windows[:,1].max()),positive_second_windows=int(np.sum(windows[:,1]>0)),windows=len(windows),
            flow_identity_max_error=a['flow_identity_max_error'],divergence_residual=float(residual)))
    a=json.loads((out/'quotient-300000/analysis.json').read_text());w=np.array([r['spectrum'] for r in a['windows']]);t=np.array([r['end'] for r in a['windows']])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.5),layout='constrained')
    fig.suptitle('Three-row motor: a weak second expanding direction remains unresolved',fontsize=14)
    axs[0].bar(t-12500,w[:,1],width=20500,color=['#db7647' if v<0 else '#218d9b' for v in w[:,1]])
    axs[0].axhline(0,color='#374151',linewidth=1);axs[0].axhline(a['spectrum'][1],color='#7c4d99',linestyle='--',label='Full-run mean')
    axs[0].set(xlabel='Original time tau',ylabel='Second non-flow rate',title='Twelve consecutive windows of 25,000 tau');axs[0].legend(frameon=False)
    axs[1].plot(t,np.cumsum(w[:,1])/np.arange(1,13),color='#7c4d99',marker='o')
    axs[1].axhline(0,color='#374151',linewidth=1);axs[1].set(xlabel='Observed time tau',ylabel='Accumulated second non-flow rate',title='Mean at 300,000 tau: +0.0028324')
    for ax in axs:ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0));ax.grid(axis='y',alpha=.15)
    fig.savefig(out/'motor-transverse-spectrum.png',dpi=180);plt.close(fig)
    (out/'analysis.json').write_text(json.dumps(dict(runs=results,
        interpretation='The autonomous flow direction is removed explicitly. The second quotient mean is positive but one of twelve long windows is negative. These finite estimates do not establish an asymptotic second positive exponent or hyperchaotic attraction.',
        independence='The 300000-unit run continues the 60000-unit endpoint after another 1000-unit warmup; these are not independent initial-condition trials.'),indent=2)+'\n')
    (out/'archive_motor_transverse.py').write_bytes(Path(__file__).read_bytes())
    (out/'README.txt').write_text('Frozen three-row motor case2742 follow-up.\nFull8D preliminary spectrum and both7Dflow-normal quotient runs are retained byte-for-byte.\nThe quotient removes the exact autonomous neutral direction, not a synchronization or cluster subspace.\nAll full-state histories are included. Time sampling every0.5tau is too coarse for a dense curve portrait.\nSource equation and priority audits remain those linked in docs/motor-asymmetry-geometry.md.\nNo new model, established hyperchaotic attractor, or historical uniqueness is claimed.\n')
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    (out/'manifest.json').write_text(json.dumps(dict(files=entries,count=len(entries),bytes=sum(v['bytes'] for v in entries)),indent=2)+'\n')
    print(json.dumps(dict(files=len(entries),results=results)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

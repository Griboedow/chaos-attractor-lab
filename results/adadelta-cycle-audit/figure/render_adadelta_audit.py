"""Visualize the stored AdaDelta cycle replay without rerunning its dynamics."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    with np.load(args.replay/'histories.npz') as d:times=d['times'];states=d['snapshots']
    with np.load(args.independent/'actual-torch.npz') as d:actual=d['states']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11})
    fig,axes=plt.subplots(1,3,figsize=(17,5.5),layout='constrained')
    steps=np.arange(1,len(actual)+1)
    axes[0].plot(steps[-500:],actual[-500:,0],color='#25738f',lw=.65)
    axes[0].set(title='Short view: apparently constant oscillation',xlabel='Actual PyTorch update',ylabel='Position x')
    axes[1].plot(steps,1-abs(actual[:,0]),color='#b95f3c',lw=2)
    axes[1].set(title='The amplitude is already changing',xlabel='Actual PyTorch update',ylabel='Amplitude loss: 1 - abs(x)')
    axes[1].ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
    axes[2].plot(times/1e6,abs(states[0,:,0]),color='#b95f3c',lw=2,label='alpha=0.1, E=1e-8')
    axes[2].plot(times/1e6,abs(states[3,:,0]),color='#25738f',lw=2,label='alpha=3, E=1e-8')
    axes[2].plot(times/1e6,abs(states[1,:,0]),color='#525b67',lw=1.5,ls='--',label='Exact controls: E=0 or alpha=2')
    assert np.array_equal(abs(states[1,:,0]),abs(states[2,:,0]))
    axes[2].set(title='Long view: substantial finite drift',xlabel='Updates (millions)',ylabel='Absolute sampled position')
    axes[2].legend(fontsize=9,loc='lower left')
    for ax in axes:ax.grid(alpha=.17)
    fig.suptitle('AdaDelta: retaining positive epsilon changes the claimed cycle',fontsize=19)
    fig.supxlabel('Source convention: curvature=1, beta=0.9, E=epsilon squared. Initial x=1, v=1; update memory d=4/alpha squared.\nThe exact positive-E period-two condition is alpha times curvature = 2. Finite drift does not prove global convergence or divergence.',fontsize=11)
    fig.savefig(out/'adadelta-cycle-drift.png',dpi=180);plt.close(fig)
    record=dict(inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.replay/'histories.npz',args.independent/'actual-torch.npz']},
        selection='Short oscillation and amplitude-loss panels use all5000 actual-library states (last500 in first panel). Long panel uses all declared saved times, with connecting guide segments. Larger-E control remains in the data but is omitted from this composition.',
        scope='Finite replay illustrated alongside the separate algebraic classification.')
    (out/'analysis.json').write_text(json.dumps(record,indent=2)+'\n');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['replay','independent','output']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())

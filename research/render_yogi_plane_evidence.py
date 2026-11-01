"""Plot stored grid outcomes at their actual coordinates and long-run evidence."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    with np.load(args.source/'histories.npz') as z:d={k:z[k] for k in z.files}
    protocol=json.loads((args.source/'protocol.json').read_text());a=json.loads((args.source/'analysis.json').read_text())
    flags=np.zeros(len(d['parameters']),dtype=int);flags[np.array(a['periods'])==0]=1;flags[a['leading_ids']]=2;flags[a['two_ids']]=3
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    for s,design in enumerate(protocol['slices']):
        ids=np.flatnonzero(d['slice_ids']==s);bs=np.array(design['beta1']);cs=np.array(design['beta2'])
        for j in range(3):
            ax=axes[s,j];z=(d['windows'][ids,:,j].mean(1) if j<2 else flags[ids]).reshape((len(cs),len(bs)))
            if j==0:
                im=ax.pcolormesh(bs,cs,z,shading='nearest',cmap='magma',vmin=0,vmax=.45)
                fig.colorbar(im,ax=ax,label='Mean first branch rate',extend='both')
            elif j==1:
                im=ax.pcolormesh(bs,cs,z,shading='nearest',cmap='RdBu_r',vmin=-.02,vmax=.02)
                fig.colorbar(im,ax=ax,label='Mean second branch rate',extend='both')
            else:
                cmap=ListedColormap(['#c9cdd3','#7868ac','#dd9550','#138f91'])
                im=ax.pcolormesh(bs,cs,z,shading='nearest',cmap=cmap,norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4))
                cb=fig.colorbar(im,ax=ax,ticks=[0,1,2,3]);cb.ax.set_yticklabels(['Detected period','Unresolved','Leading rule','Two-direction rule'])
            ax.set(xlabel='beta1',ylabel='beta2',title=f"{design['name']} / e={design['epsilon_normalized']:.3g}",
                   xlim=(max(0,bs[0]-(bs[1]-bs[0])/2),min(1,bs[-1]+(bs[1]-bs[0])/2)),
                   ylim=(max(0,cs[0]-(cs[1]-cs[0])/2),min(1,cs[-1]+(cs[1]-cs[0])/2)))
            line=np.linspace(max(bs[0],np.sqrt(cs[0])),min(bs[-1],np.sqrt(cs[-1])),100)
            ax.plot(line,line*line,color='#1d5d30',ls=':',lw=1.3)
    fig.suptitle('Yogi: two cold-start parameter slices / 2,048 complete histories',fontsize=17)
    fig.supxlabel('100,000 warmup + 200,000 observed steps per sampled setting. Dotted curve: beta1 squared = beta2.\nCells are centered on actual sampled coordinates; colors saturate at stated limits. No open chaotic region is proved.',fontsize=10)
    fig.savefig(out/'yogi-parameter-plane.png',dpi=165);fig.savefig(out/'yogi-parameter-plane.pdf');plt.close(fig)
    with np.load(args.follow/'histories.npz') as z:f={k:z[k] for k in z.files}
    fig,axes=plt.subplots(2,2,figsize=(15,9),layout='constrained')
    for row,ident in enumerate([670,125]):
        k=list(f['ids']).index(ident);w=f['windows'][k];x=(np.arange(len(w))+.5)*.05
        for j in range(2):
            ax=axes[row,j];ax.axhline(0,color='#515b6a',lw=1)
            ax.plot(x,w[:,j],color=['#cc6543','#138f91'][j],lw=1.6)
            ax.axhline(w[:,j].mean(),color='#26384a',ls='--',lw=1)
            ax.set(title=f'G{ident:04d} / direction {j+1}',xlabel='Observed steps (millions)',ylabel='Branch growth rate per step')
            ax.text(.025,.96,f'mean {w[:,j].mean():+.6f}\nminimum {w[:,j].min():+.6f}',transform=ax.transAxes,va='top',bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
            lo=min(0,w[:,j].min());hi=max(0,w[:,j].max());span=max(hi-lo,.005);ax.set_ylim(lo-.12*span,hi+.3*span)
            ax.grid(alpha=.18)
    fig.suptitle('A feathered arc and a stronger two-direction candidate',fontsize=18)
    fig.supxlabel('Fresh cold starts / 100,000 warmup + 5,000,000 observed steps / 100 disjoint windows of 50,000 steps.\nFull three-dimensional derivative on each smooth branch. Discontinuity and finite observation limit the conclusion.',fontsize=11)
    fig.savefig(out/'yogi-grid-directions.png',dpi=180);fig.savefig(out/'yogi-grid-directions.pdf');plt.close(fig)
    record=dict(inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.source/'histories.npz',args.source/'analysis.json',args.source/'protocol.json',args.follow/'histories.npz']},
        correction='Grid cells centered on actual parameter coordinates with pcolormesh, replacing the first gallery plot whose imshow extent placed outer centers at bin edges. No scientific data changed.',
        selected_evidence_ids=[670,125],scope='Finite numerical classification and windowed evidence, not a proof of hyperchaos or novelty.')
    (out/'analysis.json').write_text(json.dumps(record,indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['source','follow','output']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())

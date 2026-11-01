"""Render the complete declared finite ensemble, without PPC classification."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    a=json.loads((args.source/'analysis.json').read_text());t=np.array(a['times'])
    bg='#071521';fg='#e6eef4';muted='#9ab0bf';colors=['#52e0d4','#f6b353','#b6a0f5']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':fg,'axes.labelcolor':muted,
        'xtick.color':muted,'ytick.color':muted,'axes.edgecolor':'#345064','axes.facecolor':bg,'figure.facecolor':bg})
    fig,axes=plt.subplots(2,2,figsize=(15,9),dpi=180)
    fig.subplots_adjust(left=.075,right=.965,bottom=.19,top=.77,hspace=.34,wspace=.2)
    fig.text(.075,.94,'HOW QUICKLY DO NEARBY MOTOR TRAJECTORIES SEPARATE?',fontsize=20,weight='bold')
    fig.text(.075,.885,'Same equations, two regimes. Full-state distances; 32 saved phases at each perturbation size.',fontsize=13,color=muted)
    summaries=[]
    for si,title in enumerate(['Thin braid  |  b = 1.500','Larger regime  |  b = 1.545']):
        axes[0,si].set_title(title,pad=14)
        for j,r in enumerate([r for r in a['records'] if r['source_index']==si]):
            d=np.array(r['mean_distance']);c=np.array(r['cross_correlation']);se=np.array(r['correlation_standard_error'])
            label=f"initial distance / extent = {r['delta']:.0e}"
            axes[0,si].plot(t[1:],d[1:],color=colors[j],lw=1.3,label=label)
            axes[1,si].plot(t[1:],c[1:],color=colors[j],lw=1.2)
            axes[1,si].fill_between(t[1:],(c-se)[1:],(c+se)[1:],color=colors[j],alpha=.09,lw=0)
            summaries.append(dict(source_index=si,delta=r['delta'],final_correlation=float(c[-1]),final_mean_distance=float(d[-1]),
                late_time_mean_correlation=float(c[t>=5000].mean()),late_time_mean_definition='Arithmetic mean of saved equally spaced checkpoints5000..10000; not an independent-sample confidence interval.'))
        axes[0,si].set_yscale('log');axes[0,si].set_ylim(1e-9,2)
        axes[0,si].legend(loc='lower right',fontsize=8,facecolor=bg,edgecolor='#345064',labelcolor=fg)
        axes[1,si].axhline(0,color=muted,lw=.6);axes[1,si].set_ylim(-.35,1.45)
        axes[1,si].set_xlabel('elapsed original time units')
        for ax in axes[:,si]:ax.set_xscale('log');ax.set_xlim(1,10000);ax.grid(alpha=.09)
    axes[0,0].set_ylabel('mean distance / source RMS extent')
    axes[1,0].set_ylabel('centered cross-correlation')
    fig.text(.075,.102,'Shading: one standard error across these 32 phases; the phases are not independent natural-measure draws.',fontsize=10,color=muted)
    fig.text(.075,.072,'The braid loses correlation more slowly. This finite test does not establish partially predictable chaos or a new mechanism.',fontsize=10,color=muted)
    fig.savefig(out/'motor-pair-dispersion.png');plt.close(fig)
    (out/'render.json').write_text(json.dumps(dict(source=str(args.source),analysis_sha256=hashlib.sha256((args.source/'analysis.json').read_bytes()).hexdigest(),
        summaries=summaries,curves='All declared checkpoints; no smoothing, fit, discarded phase or automatic classification.',prior='https://doi.org/10.1038/s41598-017-01083-x'),indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

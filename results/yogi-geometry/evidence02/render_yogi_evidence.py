"""Compose the observed geometry and independent expansion checks for Yogi4214."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);root=args.external
    source=root/'yogi-spectra-five-million/histories.npz';d=np.load(source);i=list(d['ids']).index(4214);w=d['windows'][i]
    upstream=json.loads((root/'yogi-upstream-analysis02/analysis.json').read_text())
    u=next(r for r in upstream['runs'] if r['id']==4214)
    precision=json.loads((root/'yogi-precision/analysis.json').read_text())
    pp=[r for r in precision['runs'] if r['id']==4214]
    assert w.shape==(100,3) and np.min(w[:,:2])>0
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(15,9),facecolor='#f4f6f8')
    gs=fig.add_gridspec(2,2,width_ratios=[1.05,1],height_ratios=[1,1],left=.035,right=.975,bottom=.16,top=.84,wspace=.23,hspace=.5)
    ax=fig.add_subplot(gs[:,0]);im=Image.open(root/'yogi-previews/case-4214-view-1.png');ax.imshow(im);ax.axis('off')
    ax.set_title('Full normalized state (x, m, v)\n2,000,000 actual iterates; rigid orthographic view',fontsize=12,pad=15)
    a=fig.add_subplot(gs[0,1]);t=np.arange(1,101)*.05
    a.plot(t,w[:,0],color='#9d5a18',lw=1,label='First direction')
    a.plot(t,w[:,1],color='#176b87',lw=1,label='Second direction')
    a.axhline(0,color='#444',lw=.8);a.set(xlabel='Observed iterations (millions)',ylabel='Branch growth per step',ylim=(-.003,.058),title='100 consecutive windows of 50,000 steps')
    a.legend(frameon=False,ncol=2,fontsize=10);a.grid(alpha=.15)
    b=fig.add_subplot(gs[1,1]);labels=['5 million\nbinary64','70 digits\nfinite separation','90 digits\nfinite separation','Actual Optax\n190,000 steps']
    windows=[w[:,1]]+[np.array(r['windows'])[:,1] for r in pp]+[np.array(u['windows'])[:,1]]
    means=[x.mean() for x in windows];lows=[x.min() for x in windows];highs=[x.max() for x in windows]
    b.errorbar(np.arange(4),means,yerr=[np.array(means)-lows,np.array(highs)-means],fmt='o',color='#176b87',capsize=6,lw=2)
    b.axhline(0,color='#444',lw=.8);b.set(xticks=np.arange(4),xticklabels=labels,ylabel='Second growth rate',title='Means and full observed window ranges',xlim=(-.5,3.5));b.grid(axis='y',alpha=.15)
    b.tick_params(axis='x',labelsize=9)
    fig.suptitle('Yogi4214: numerical evidence for two expanding directions',fontsize=19)
    fig.text(.5,.035,'Window lengths differ: 50,000 (binary64), 10,000 (70/90 digits and Optax). Ranges are not confidence intervals.\nOne short Optax window is negative. Finite observations do not prove asymptotic hyperchaotic attraction.',ha='center',fontsize=10,color='#3c4859')
    fig.savefig(out/'yogi-two-directions.png',dpi=180);plt.close(fig)
    record=dict(id=4214,source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),labels=labels,
        second_means=means,second_window_minima=lows,second_window_maxima=highs,
        image_source=str(root/'yogi-previews/case-4214-view-1.png'),
        image_source_sha256=hashlib.sha256((root/'yogi-previews/case-4214-view-1.png').read_bytes()).hexdigest(),
        scope='No new simulation or geometric transformation; observed window ranges, not error bars on an asymptotic estimate.')
    (out/'analysis.json').write_text(json.dumps(record,indent=2)+'\n');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

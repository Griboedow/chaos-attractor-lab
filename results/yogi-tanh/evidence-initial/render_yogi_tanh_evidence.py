"""Compare finite spectra across long runs, independent precision and actual library states."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    def read(name):return json.loads((args.external/name/'analysis.json').read_text())
    long=read('yogi-tanh-five-million');precision=read('yogi-tanh-precision');cohorts=read('yogi-tanh-cohorts')
    upstream=read('yogi-tanh-upstream-analysis02')['runs']+read('yogi-tanh-upstream-strong-analysis')['runs']
    with np.load(args.external/'yogi-tanh-five-million/histories.npz') as f:d={k:f[k] for k in f.files}
    colors=['#ffbd69','#68d5eb'];bg='#071320';fg='#e8f0f5';muted='#a4b9ca'
    fig,axes=plt.subplots(2,3,figsize=(18,8.5),facecolor=bg,gridspec_kw={'width_ratios':[1.25,1.1,1.05]});records=[]
    for row,ident in enumerate([39,7239]):
        i=list(d['ids']).index(ident);w=d['windows'][i];x=np.arange(1,len(w)+1)*50000/1e6
        for ax in axes[row]:
            ax.set_facecolor(bg);ax.tick_params(colors=muted,labelsize=9)
            for spine in ax.spines.values():spine.set_color('#294153')
        ax=axes[row,0]
        for j in [0,1]:ax.plot(x,w[:,j],color=colors[j],lw=1.2,label=f'growth rate {j+1}')
        ax.axhline(0,color=muted,lw=.8);ax.set_ylim(bottom=-.006);ax.set_title(f'S{ident:04d} / five million steps',color=fg,fontsize=14)
        ax.set_xlabel('Observed steps (millions); windows of 50,000',color=muted,fontsize=9);ax.set_ylabel('Natural log growth / step',color=muted,fontsize=9)
        ax.legend(frameon=False,labelcolor=fg,fontsize=9,loc='upper right')
        datasets=[w[:,1]]
        for digits in [70,90]:datasets.append(np.array(next(r for r in precision['runs'] if r['id']==ident and r['digits']==digits)['windows'])[:,1])
        original=next(r for r in upstream if r['id']==ident)
        folder='yogi-tanh-upstream-analysis02' if ident==39 else 'yogi-tanh-upstream-strong-analysis'
        with np.load(args.external/folder/'windows.npz') as f:datasets.append(f['windows'][list(f['ids']).index(ident),:,1])
        ax=axes[row,1]
        for col,y in enumerate(datasets):
            ax.scatter(np.full(len(y),col),y,s=13,color=colors[1],alpha=.5,zorder=3)
            ax.plot([col-.16,col+.16],[y.mean(),y.mean()],color=colors[0],lw=3,zorder=4)
            records.append(dict(id=ident,method=['binary64_5m','independent_70digits','independent_90digits','actual_tensorflow'][col],windows=y.tolist(),minimum=float(y.min()),mean=float(y.mean())))
        ax.axhline(0,color=muted,lw=.8);ax.set_ylim(bottom=-.0015);ax.set_xticks(range(4),['5M steps','70 digits','90 digits','TensorFlow']);ax.set_title('Second growth rate / every window',color=fg,fontsize=13)
        ax.set_xlabel('Dots = windows; gold line = mean',color=muted,fontsize=9)
        ax=axes[row,2];kinds=['autonomous_cold','library_style_cold','parameter_neighbor','state_box'];labels=['Cold starts','Bias starts','Nearby parameters','State box']
        selected=[next(r for r in cohorts['summary'] if r['id']==ident and r['kind']==kind) for kind in kinds]
        for j,r in enumerate(selected):
            ax.barh(j,1,color='#223849',height=.55);ax.barh(j,r['two']/r['total'],color=colors[1],height=.55)
            ax.text(.5,j,f"{r['two']} / {r['total']}",ha='center',va='center',color=bg,weight='bold',fontsize=11)
        ax.set_yticks(range(4),labels);ax.invert_yaxis();ax.set_xlim(0,1.04);ax.set_xticks([0,.5,1],['0','1/2','all']);ax.set_title('Two positive rates / fresh runs',color=fg,fontsize=13)
        ax.set_xlabel('Declared runs passing all window thresholds',color=muted,fontsize=9)
    fig.text(.045,.96,'SMOOTH YOGI / TWO EXPANDING DIRECTIONS',fontsize=23,color=fg,weight='bold')
    fig.text(.045,.91,'The tanh option already exists in TensorFlow Addons. These are finite numerical checks on a scalar quadratic.',fontsize=12,color=muted)
    fig.text(.045,.055,'Independent precision uses finite separations, without the analytic derivative. TensorFlow uses unmodified TFA 0.22.0 at its recorded effective coefficients.',fontsize=10,color=muted)
    fig.text(.045,.025,'Precision and TensorFlow windows: 10,000 steps. Fresh runs: 100,000 warm-up + 500,000 observed. Counts are evidence of attraction, not basin probabilities or a chaos theorem.',fontsize=10,color=muted)
    fig.subplots_adjust(left=.055,right=.975,bottom=.14,top=.83,wspace=.39,hspace=.42)
    fig.savefig(out/'yogi-tanh-evidence.png',dpi=180,facecolor=bg);plt.close(fig)
    (out/'analysis.json').write_text(json.dumps(dict(comparisons=records,source_directories=['yogi-tanh-five-million','yogi-tanh-precision','yogi-tanh-cohorts','yogi-tanh-upstream-analysis02','yogi-tanh-upstream-strong-analysis']),indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

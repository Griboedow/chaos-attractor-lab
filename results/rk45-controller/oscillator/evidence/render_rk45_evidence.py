"""Evidence panels and initial-volume approach, directly from retained arrays."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from render_twist import camera

COLORS=['#007d92','#a64f00','#8c45aa']
IDS=[335,448,2006]

def read(path):
    with np.load(path) as f:return {k:f[k] for k in f.files}

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);base=args.external
    long=read(base/'rk45-oscillator-long01/follow.npz');scout=read(base/'rk45-oscillator-scout01/screen.npz')
    cohort=read(base/'rk45-oscillator-cohort01/cohort.npz');records=[]
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(15,10))
    ax=axes[0,0]
    for i,(ident,color) in enumerate(zip(IDS,COLORS)):
        r=long['rates'][list(long['ids']).index(ident)];ax.plot(np.arange(1,21)*.25,r[:,0],color=color,marker='o',ms=3,label=f'D{ident}')
    ax.axhline(0,color='#76818c',lw=.8);ax.set(xlabel='Observed accepted steps (millions)',ylabel='Leading finite growth rate / step',title='A  Five million steps from a fresh start');ax.legend(ncol=3,frameon=False)
    ax=axes[0,1]
    methods=['70 digits','90 digits','Actual SciPy'];xx=np.arange(3)
    for i,(ident,color) in enumerate(zip(IDS,COLORS)):
        means=[];low=[];high=[]
        for digits in [70,90]:
            d=json.loads((base/f'rk45-oscillator-precision01/case-{ident}-{digits}.json').read_text())
            # Saved independent finite-separation estimates in actual q,p,h.
            key='rates' if 'rates' in d else 'windows'
            rr=np.array(d[key]);means.append(rr[:,0].mean());low.append(rr[:,0].min());high.append(rr[:,0].max())
        d=read(base/f'rk45-oscillator-upstream01/case-{ident}-analysis.npz')['rates'];means.append(d[:,0].mean());low.append(d[:,0].min());high.append(d[:,0].max())
        means=np.array(means);ax.errorbar(xx+(i-1)*.16,means,yerr=[means-low,high-means],fmt='o',capsize=4,color=color,label=f'D{ident}')
        records.append(dict(id=ident,independent_methods=methods,means=means.tolist(),minima=low,maxima=high))
    ax.set(xticks=xx,xticklabels=methods,ylabel='Leading finite growth rate / step',title='B  Independent checks');ax.axhline(0,color='#76818c',lw=.8)
    ax.text(.02,.23,'Bars span all measured windows, not confidence intervals.',transform=ax.transAxes,fontsize=10,color='#5a6470')
    ax=axes[1,0]
    bad=~scout['flags'];p=scout['params'];ax.scatter(p[bad,0],p[bad,1],s=9,c='#b7c2c8',alpha=.65,label='Not flagged')
    im=ax.scatter(p[~bad,0],p[~bad,1],s=13,c=scout['rates'][~bad].mean(1),cmap='viridis',vmin=0,vmax=.4)
    for ident,color in zip(IDS,COLORS):ax.scatter(*p[ident],s=110,facecolors='none',edgecolors=color,lw=1.7)
    ax.set(yscale='log',xlabel='Damping',ylabel='Relative tolerance',title='C  All 2,051 screened settings');fig.colorbar(im,ax=ax,label='Scout leading growth rate',pad=.025)
    ax=axes[1,1]
    groups=['start','step','neighbor'];labels=['24 solution starts','6 initial steps','32 nearby parameters'];x=np.arange(3)
    for i,(ident,color) in enumerate(zip(IDS,COLORS)):
        passed=[];total=[]
        for group in groups:
            mask=(cohort['ids']==ident)&(cohort['groups']==group);passed.append(int(np.sum(cohort['flags'][mask])));total.append(int(np.sum(mask)))
        bars=ax.bar(x+(i-1)*.24,np.array(passed)/total,width=.22,color=color,label=f'D{ident}')
        for bar,n,t in zip(bars,passed,total):ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+.02,f'{n}/{t}',ha='center',fontsize=9)
    ax.set(xticks=x,xticklabels=labels,ylim=(0,1.23),ylabel='Fraction passing the declared finite criterion',title='D  Nearby periodic outcomes are retained');ax.legend(ncol=3,loc='upper left',frameon=False)
    fig.suptitle('RK45 oscillator candidates / reproducibility and limits',fontsize=20,x=.055,ha='left')
    fig.text(.055,.025,'The ODE is stable. Rates describe its numerical accepted-step map; finite positive estimates do not certify chaos or historical novelty.',fontsize=11,color='#46515e')
    fig.subplots_adjust(top=.9,bottom=.11,left=.07,right=.97,hspace=.37,wspace=.31)
    fig.savefig(args.out/'rk45-oscillator-evidence.png',dpi=180);fig.savefig(args.out/'rk45-oscillator-evidence.pdf');plt.close(fig)
    # Equal-scale projections of unmodified state coordinates in each row.
    render=json.loads((base/'rk45-oscillator-render01/analysis.json').read_text())['runs']
    fig,axes=plt.subplots(2,4,figsize=(18,9),gridspec_kw={'width_ratios':[1,1,1,1.25]})
    approach=[]
    for row,ident in enumerate([335,448]):
        d=read(base/f'rk45-oscillator-approach01/case-{ident}.npz');r=next(x for x in render if x['id']==ident)
        basis=np.array(r['basis']);center=np.array(r['center']);proj=camera(r['azimuth'],r['elevation'])
        coords=(d['snapshots']-center)@basis@proj.T;half=float(np.nanmax(np.abs(coords)))*1.03
        for col,t in enumerate([0,16,1024]):
            k=list(d['times']).index(t);ax=axes[row,col];xy=coords[k];good=np.all(np.isfinite(xy),axis=1)
            ax.scatter(xy[good,0],xy[good,1],s=.25,c=COLORS[row*2],alpha=.35,rasterized=True)
            ax.set(xlim=(-half,half),ylim=(-half,half),aspect='equal',title=f'D{ident} / step {t:,}');ax.set_axis_off()
        ax=axes[row,3];quant=np.nanquantile(d['distances'],[.1,.5,.9],axis=1);times=d['times'];ax.plot(times,quant[1],'-o',color=COLORS[row*2],ms=4);ax.fill_between(times,quant[0],quant[2],color=COLORS[row*2],alpha=.18)
        floor=float(np.median(d['reference_split_distances']));ax.axhline(floor,color='#777777',ls='--',label='Reference sampling scale')
        ax.set(xscale='symlog',yscale='log',xlim=(0,1150),xlabel='Accepted steps',ylabel='Distance to reference cloud',title=f'D{ident} / approach distance');ax.set_xticks([0,4,16,64,256,1024],['0','4','16','64','256','1024']);ax.legend(fontsize=8,frameon=False)
        approach.append(dict(id=ident,common_halfwidth=half,projected_count_per_snapshot=np.isfinite(coords).all(2).sum(1).tolist(),reference_sampling_median=floor))
    fig.suptitle('32,768 initial states approach each selected cloud',fontsize=21,x=.03,ha='left')
    fig.text(.03,.035,'Each row uses one rigid camera and one scale. Box: q/atol and p/atol in [-4,4], h in [0.25,8]. All initial points retained.',fontsize=11)
    fig.text(.03,.008,'Distance bands show the 10th–90th percentiles to 200,000 reference points; the dashed line is a finite sampling benchmark, not zero.',fontsize=11,color='#46515e')
    fig.subplots_adjust(top=.88,bottom=.12,left=.02,right=.98,hspace=.32,wspace=.24)
    fig.savefig(args.out/'rk45-oscillator-approach.png',dpi=180);fig.savefig(args.out/'rk45-oscillator-approach.pdf');plt.close(fig)
    (args.out/'report.json').write_text(json.dumps(dict(evidence=records,approach=approach,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2),encoding='utf-8')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print('Rendered evidence and approach',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--external',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

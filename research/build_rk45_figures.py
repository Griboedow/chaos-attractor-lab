"""Publication figures using only the frozen RK45 arrays."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify_rk23_controller import read,sha
from render_twist import colorize,camera

COLORS=['#007d92','#ad570d','#8c45aa']


def load(root,relative,sources):
    path=root/relative;sources.append(dict(path='results/rk45-controller/'+relative,sha256=sha(path)))
    return json.loads(path.read_text()) if path.suffix=='.json' else read(path)


def finish(fig,out,name):
    fig.savefig(out/(name+'.png'),dpi=300,facecolor='white');fig.savefig(out/(name+'.pdf'),facecolor='white');plt.close(fig)


def figures(root,out):
    out.mkdir(parents=True,exist_ok=True);sources=[]
    plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    hs=load(root,'oscillator/render/portrait-histograms.npz',sources)['histograms'];render=load(root,'oscillator/render/analysis.json',sources)['runs']
    fig,axes=plt.subplots(1,2,figsize=(7,2.5))
    for ax,h,row in zip(axes,hs,render):
        im,_=colorize(h,row['log_density_vmax']);ax.imshow(im,interpolation='nearest',origin='upper');ax.set_axis_off();ax.set_title(row['title'].replace(' / ',' | '),pad=5)
    fig.subplots_adjust(left=.005,right=.995,bottom=.01,top=.88,wspace=.035);finish(fig,out,'figure-1-geometry')
    long=load(root,'oscillator/five-million/follow.npz',sources);rates=long['rates'][list(long['ids']).index(335)]
    fig,axes=plt.subplots(1,2,figsize=(7,2.7));ax=axes[0]
    for j,color in enumerate(COLORS[:2]):ax.plot(np.arange(1,21)*.25,rates[:,j],'-o',ms=2,lw=.8,color=color,label=f'Direction {j+1}')
    ax.axhline(0,color='#7e8c95',lw=.6);ax.set(xlabel='Observed steps (millions)',ylabel='Finite growth / step',title='D335: two expanding directions');ax.legend(fontsize=7,frameon=False,loc='best')
    ax=axes[1];series=[rates]
    for digits in [70,90]:series.append(np.array(load(root,f'oscillator/precision/case-335-{digits}.json',sources)['rates']))
    series.append(load(root,'oscillator/upstream/case-335-analysis.npz',sources)['rates'])
    for j,color in enumerate(COLORS[:2]):
        means=np.array([v[:,j].mean() for v in series]);lower=np.array([v[:,j].min() for v in series]);upper=np.array([v[:,j].max() for v in series])
        ax.errorbar(np.arange(4)+(j-.5)*.15,means,yerr=[means-lower,upper-means],fmt='o',ms=3,capsize=2,lw=.8,color=color)
    ax.axhline(0,color='#7e8c95',lw=.6);ax.set(xticks=np.arange(4),xticklabels=['5M','70 digits','90 digits','SciPy'],title='Independent checks',ylabel='Finite growth / step')
    fig.subplots_adjust(left=.1,right=.99,bottom=.19,top=.86,wspace=.43);finish(fig,out,'figure-2-evidence')
    fig,axes=plt.subplots(2,3,figsize=(7,3.7))
    for row,ident in enumerate([335,448]):
        d=load(root,f'oscillator/approach/case-{ident}.npz',sources);r=next(v for v in render if v['id']==ident)
        coords=(d['snapshots']-np.array(r['center']))@np.array(r['basis'])@camera(r['azimuth'],r['elevation']).T;half=np.max(np.abs(coords))*1.03
        for col,t in enumerate([0,16,1024]):
            xy=coords[list(d['times']).index(t)];ax=axes[row,col];ax.scatter(xy[:,0],xy[:,1],s=.12,c=COLORS[row*2],alpha=.35,rasterized=True);ax.set(xlim=(-half,half),ylim=(-half,half),aspect='equal',title=f'D{ident} / step {t:,}');ax.set_axis_off()
    fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.92,wspace=.02,hspace=.2);finish(fig,out,'figure-3-approach')
    d=load(root,'oscillator/damping/scan.npz',sources);base=load(root,'oscillator/damping/report.json',sources)['base_params'][0]
    fig,axes=plt.subplots(2,2,figsize=(7,4.7))
    for row,group in enumerate(['wide','narrow']):
        m=d['groups']==group;x=d['params'][m,0];tail=d['tails'][m,:,1];ax=axes[row,0]
        ax.scatter(np.repeat(x,tail.shape[1]),tail.reshape(-1),s=.015,c=COLORS[0],alpha=.2,rasterized=True)
        ax.axvline(base,color='#ad570d',lw=.7,ls='--');ax.set(xlabel='Damping',ylabel='Late velocity / atol',title=f'{group.capitalize()} slice: every retained tail point')
        ax=axes[row,1]
        for j,color in enumerate(COLORS[:2]):
            rr=d['rates'][m,:,j];ax.plot(x,rr.mean(1),color=color,lw=.7,label=f'Direction {j+1}');ax.fill_between(x,rr.min(1),rr.max(1),color=color,alpha=.2,lw=0)
        ax.axhline(0,color='#7e8c95',lw=.6);ax.axvline(base,color='#777777',lw=.7,ls='--');ax.set(xlabel='Damping',ylabel='Finite growth / step',title=f'{group.capitalize()} slice: all window ranges')
        if row==0:ax.legend(fontsize=6,frameon=False,loc='best')
        for ax in axes[row]:ax.ticklabel_format(axis='x',useOffset=False,style='plain')
    fig.subplots_adjust(left=.1,right=.985,bottom=.11,top=.94,wspace=.39,hspace=.53);finish(fig,out,'figure-4-damping')
    (out/'provenance.json').write_text(json.dumps(dict(sources=sources,scope='Publication layouts from frozen exact histograms, finite-window arrays and volume snapshots; no new simulation or geometric deformation.'),indent=2),encoding='utf-8')

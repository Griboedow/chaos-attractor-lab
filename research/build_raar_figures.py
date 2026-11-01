"""Publication figures from the frozen spatial-circle packet."""
import json,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify_rk23_controller import read,sha
from render_twist import colorize,camera


def figures(root,out):
    out.mkdir(parents=True,exist_ok=True);sources=[]
    def load(relative):
        path=root/relative;sources.append(dict(path='results/raar-circles/'+relative,sha256=sha(path)))
        return json.loads(path.read_text()) if path.suffix=='.json' else read(path)
    def finish(fig,name):
        fig.savefig(out/(name+'.png'),dpi=300,facecolor='white');fig.savefig(out/(name+'.pdf'),facecolor='white');plt.close(fig)
    plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False})
    h=load('render/portrait-histograms.npz')['histograms'];render=load('render/analysis.json')['runs']
    fig,axes=plt.subplots(1,2,figsize=(7,2.5))
    for ax,hh,r in zip(axes,h,render):
        im,_=colorize(hh,r['log_density_vmax']);ax.imshow(im,interpolation='nearest');ax.set_axis_off();ax.set_title(r['title'].replace(' / ',' | '))
    fig.subplots_adjust(left=.005,right=.995,bottom=.01,top=.88,wspace=.035);finish(fig,'figure-1-geometry')
    long=load('five-million/follow.npz');fig,axes=plt.subplots(1,2,figsize=(7,2.7))
    for i,col in [(1051,'#007d92'),(1483,'#ad570d'),(820,'#8c45aa')]:
        r=long['rates'][list(long['ids']).index(i)];axes[0].plot(np.arange(1,21)*.25,r[:,0],'-o',ms=2,lw=.8,c=col,label=f'R{i}')
    axes[0].set(xlabel='Observed steps (millions)',ylabel='Leading finite growth / step',title='Five-million-step checks');axes[0].legend(ncol=3,fontsize=6,frameon=False)
    axes[0].axhline(0,c='#777777',lw=.6)
    series=[long['rates'][list(long['ids']).index(820)]]
    for digits in [70,90]:series.append(np.array(load(f'precision/case-820-{digits}.json')['rates']))
    series.append(load('upstream/case-820.npz')['rates'])
    m=np.array([r[:,1].mean() for r in series]);lo=np.array([r[:,1].min() for r in series]);hi=np.array([r[:,1].max() for r in series])
    axes[1].errorbar(np.arange(4),m,yerr=[m-lo,hi-m],fmt='o',ms=4,capsize=3,color='#8c45aa');axes[1].axhline(0,c='#777777',lw=.6)
    axes[1].set(xticks=np.arange(4),xticklabels=['5M','70 digits','90 digits','Library'],ylabel='Second finite growth / step',title='R820: positive means, one negative window')
    fig.subplots_adjust(left=.1,right=.99,bottom=.19,top=.85,wspace=.43);finish(fig,'figure-2-evidence')
    scout=load('scout/screen.npz');fig=plt.figure(figsize=(7,3.8));t=np.linspace(0,2*np.pi,600);A=np.c_[np.cos(t),np.sin(t),np.zeros_like(t)]
    for j,i in enumerate([1051,1483,820,1679]):
        p=scout['params'][i];b,r,theta=p[:3];B=p[3:]+r*np.c_[np.cos(t),np.cos(theta)*np.sin(t),np.sin(theta)*np.sin(t)]
        ax=fig.add_subplot(2,2,j+1,projection='3d');ax.plot(*A.T,c='#007d92',lw=1.2,label='A');ax.plot(*B.T,c='#ad570d',lw=1.2,label='B')
        cloud=np.r_[A,B];mid=(cloud.min(0)+cloud.max(0))/2;half=(cloud.max(0)-cloud.min(0)).max()*.54
        ax.set(xlim=(mid[0]-half,mid[0]+half),ylim=(mid[1]-half,mid[1]+half),zlim=(mid[2]-half,mid[2]+half),title=f'R{i}: input-circle linking number '+('−1' if i==1679 else '0'))
        ax.set_box_aspect((1,1,1));ax.view_init(24,-48);ax.set_axis_off()
        if j==0:ax.legend(loc='upper left',frameon=False,fontsize=7)
    fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.94,wspace=.03,hspace=.05);finish(fig,'figure-3-constraints')
    fig,axes=plt.subplots(2,3,figsize=(7,2.7))
    for row,i in enumerate([1051,1483]):
        d=load(f'approach/case-{i}.npz');r=next(v for v in render if v['id']==i);xy=(d['snapshots']-np.array(r['center']))@np.array(r['basis'])@camera(r['azimuth'],r['elevation']).T;half=np.max(np.abs(xy))*1.03
        for col,step in enumerate([0,16,1024]):
            z=xy[list(d['times']).index(step)];ax=axes[row,col];ax.scatter(z[:,0],z[:,1],s=.12,c=['#007d92','#8c45aa'][row],alpha=.35,rasterized=True);ax.set(xlim=(-half,half),ylim=(-half,half),aspect='equal',title=f'R{i} / step {step:,}');ax.set_axis_off()
    fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.92,wspace=.02,hspace=.2);finish(fig,'figure-4-approach')
    for origin,name in [('slice-figures/raar-relaxation','figure-5-relaxation'),('slice-figures/raar-library-clouds','figure-6-library')]:
        for suffix in ['.pdf','.png']:
            source=root/(origin+suffix);sources.append(dict(path='results/raar-circles/'+origin+suffix,sha256=sha(source)));shutil.copy2(source,out/(name+suffix))
    (out/'sources.json').write_text(json.dumps(dict(sources=sources),indent=2),encoding='utf-8')

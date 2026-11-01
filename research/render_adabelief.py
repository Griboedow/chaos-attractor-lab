"""Literal state plots of the recorded discrete map; no interpolated orbit."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.animation import FuncAnimation,PillowWriter

BG='#0b1220';FG='#edf3fb';MUTED='#9eafc4';CYAN='#55d8df';ORANGE='#ffa76d'
PALETTE=LinearSegmentedColormap.from_list('visible_magma',plt.colormaps['magma'](np.linspace(.23,1,256)))


def style(ax,three=False):
    ax.set_facecolor(BG);ax.tick_params(colors=MUTED,labelsize=8)
    for spine in ax.spines.values():spine.set_color('#40516b')
    ax.xaxis.label.set_color(MUTED);ax.yaxis.label.set_color(MUTED)
    if three:
        ax.zaxis.label.set_color(MUTED)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.pane.fill=False;axis.pane.set_edgecolor('#314057')
            axis._axinfo['grid']['color']=(.28,.35,.46,.25)
    else:ax.grid(alpha=.14,color=MUTED)


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--gif',action='store_true')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    primary=np.load(a.source/'long/run-000.npz');secondary=np.load(a.source/'long/run-001.npz')
    z=primary['tail'];w=primary['windows']
    fig=plt.figure(figsize=(14,7.5),facecolor=BG)
    gs=fig.add_gridspec(2,2,width_ratios=[1.3,1],left=.04,right=.97,bottom=.11,top=.85,wspace=.18,hspace=.34)
    ax=fig.add_subplot(gs[:,0],projection='3d')
    ax.scatter(*z.T,c=z[:,2],cmap=PALETTE,s=.11,alpha=.65,rasterized=True)
    ax.set(xlabel='x = parameter / learning rate',ylabel='m = scaled gradient memory',zlabel='s = scaled error memory')
    ax.view_init(23,-64);style(ax,True)
    ax=fig.add_subplot(gs[0,1]);ax.scatter(z[:,0],z[:,1],c=z[:,2],cmap=PALETTE,s=.09,alpha=.55,rasterized=True)
    ax.set(xlabel='x',ylabel='m',title='The same data, viewed in two coordinates');ax.title.set_color(FG);style(ax)
    ax=fig.add_subplot(gs[1,1]);tt=np.arange(1,len(w)+1)*10000
    ax.plot(tt,w[:,0],color=ORANGE,lw=.75,alpha=.65,label='10,000-step windows')
    ax.plot(tt,np.cumsum(w[:,0])/np.arange(1,len(w)+1),color=CYAN,lw=2,label='cumulative estimate')
    ax.axhline(0,color=MUTED,lw=.8);ax.set(xlabel='observed iterations',ylabel='leading exponent per iteration',ylim=(-.006,.225))
    ax.legend(facecolor=BG,labelcolor=FG,edgecolor='#314057',fontsize=9);style(ax)
    fig.text(.05,.95,'CHAOS INSIDE AN OPTIMIZER',color=CYAN,fontsize=11,weight='bold')
    fig.text(.05,.897,'A single convex bowl. Three internal states. Persistent numerical chaos.',color=FG,fontsize=18,weight='bold')
    fig.text(.05,.035,'AdaBelief on f(theta) = theta^2 / 2  |  beta1 = 0.52, beta2 = 0.15  |  2,000,000 observed steps\n200,000 unthinned states; color = s; separate linear axis scales. Finite-time evidence; novelty and asymptotic attraction unresolved.',color=MUTED,fontsize=9,linespacing=1.6)
    fig.savefig(a.output/'adabelief-candidate.png',dpi=180,facecolor=BG);plt.close(fig)
    pilot=np.load(a.source/'pilot.npz')
    fig,axes=plt.subplots(1,3,figsize=(15,4.5),facecolor=BG)
    sc=axes[0].scatter(pilot['b'],pilot['c'],c=np.clip(pilot['windows'].min(0),-.02,.22),s=6,cmap='coolwarm',vmin=-.02,vmax=.22)
    axes[0].scatter([.52,.9],[.15,.99],marker='x',s=60,c='black')
    axes[0].set(xlabel='beta1',ylabel='beta2',title='4,101 exploratory lanes')
    cb=fig.colorbar(sc,ax=axes[0]);cb.ax.tick_params(colors=MUTED,labelsize=8)
    cb.set_label('minimum window growth (clipped)',color=MUTED)
    for zz,col,label in [(z,ORANGE,'branching candidate'),(secondary['tail'],CYAN,'near-standard memory')]:
        axes[1].plot(np.arange(300),zz[-300:,0],lw=.6,color=col,label=label)
    axes[1].set(xlabel='last 300 iterations',ylabel='x',title='Deterministic discrete updates');axes[1].legend(facecolor=BG,labelcolor=FG,fontsize=8)
    eps=json.loads((a.source/'epsilon/summary.json').read_text())['runs']
    for start,col,label in [(0,ORANGE,'beta1=.52, beta2=.15'),(5,CYAN,'beta1=.9, beta2=.99'),(10,'#b78cf4','beta1=.9, beta2=.999')]:
        rows=eps[start:start+5]
        axes[2].semilogx([r['parameters'][2]*1e-4 for r in rows],[r['spectrum'][0] for r in rows],'-o',ms=4,color=col,label=label)
    axes[2].axhline(0,color=MUTED,lw=.7);axes[2].set(xlabel='original epsilon, learning rate = .01',ylabel='leading exponent',title='Regularization can remove the chaos')
    axes[2].legend(facecolor=BG,labelcolor=FG,fontsize=8)
    for ax in axes:style(ax);ax.title.set_color(FG)
    fig.tight_layout();fig.savefig(a.output/'adabelief-evidence.png',dpi=170,facecolor=BG);plt.close(fig)
    if a.gif:
        fig=plt.figure(figsize=(7.8,6),facecolor=BG);ax=fig.add_subplot(projection='3d')
        points=z[::4] # fixed stride stated in caption, no geometric filtering
        ax.scatter(*points.T,c=points[:,2],cmap=PALETTE,s=.2,alpha=.7)
        ax.set(xlabel='x',ylabel='m',zlabel='s');style(ax,True)
        fig.text(.08,.94,'Scalar AdaBelief: full optimizer state',color=FG,fontsize=15,weight='bold')
        fig.text(.08,.035,'Fixed data; rotating camera. Every fourth saved state; color = s.\nNumerical chaos candidate; no claim of a new topological class.',color=MUTED,fontsize=9)
        def update(i):ax.view_init(23,-64+i*5);return []
        animation=FuncAnimation(fig,update,frames=72,interval=90)
        animation.save(a.output/'adabelief-rotation.gif',writer=PillowWriter(fps=11),dpi=110)
        plt.close(fig)


if __name__=='__main__':main()

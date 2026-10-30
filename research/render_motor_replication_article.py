"""Vector scientific figures and a rotating view for the motor manuscript."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
COLORS=['#007f8b','#ce6738','#7754a0','#2f8b60','#c19925']
INK='#172e43'
MUTED='#596d7d'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14,'text.color':INK,
                         'axes.labelcolor':INK,'axes.edgecolor':'#abb7c0',
                         'xtick.color':MUTED,'ytick.color':MUTED,'axes.titlesize':14,
                         'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
                         'figure.facecolor':'white','axes.facecolor':'white',
                         'savefig.facecolor':'white'})


def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.15,linewidth=.6)


def panel(fig,x,y,label,title):
    fig.text(x,y,label,fontsize=18,fontweight='bold',va='center')
    fig.text(x+.035,y,title,fontsize=14,va='center')


def save(fig,out,name,paths):
    for suffix in ('pdf','svg','png'):
        path=out/f'{name}.{suffix}'
        fig.savefig(path,dpi=220)
        paths.append(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'papers/replication-instability/figures')
    parser.add_argument('--animate',action='store_true')
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    attraction=ROOT/'results/motor-cycle-attraction-visual/local-0p03'
    nonlinear=ROOT/'results/motor-replication-article/nonlinear'
    a=json.loads((attraction/'analysis.json').read_text())
    n=json.loads((nonlinear/'analysis.json').read_text())
    p=json.loads((nonlinear/'protocol.json').read_text())
    index=json.loads((ROOT/'results/motor-replication-article/exact-index-recheck.json').read_text())
    assert a['all_six_approach_same_reference_orbit'] and n['status']=='passed'
    assert index['status']=='exact_replication_index_certified'
    oldpaths=[attraction/f'trajectory-{i:02}.npz' for i in range(6)]
    newpaths=[nonlinear/f'trajectory-{i:02}.npz' for i in range(8)]
    for path,record in zip(oldpaths,a['runs']):
        assert sha(path)==record['trajectory_sha256']
    for path,record in zip(newpaths,n['runs']):
        assert sha(path)==record['trajectory_sha256']
    histories=[np.load(path) for path in oldpaths]
    new=[np.load(path) for path in newpaths]
    reference=np.load(attraction/'reference.npz')
    period=float(reference['period']);time=reference['time']/period
    cycle=reference['states'];orbit=cycle[::6,[16,17,20]]
    style();outputs=[]

    # The coordinates are actual row displacements, with no fitted embedding.
    fig=plt.figure(figsize=(10.8,5.0))
    panel(fig,.045,.946,'A','Physical orbit projection')
    ax=fig.add_axes([.015,.12,.49,.72],projection='3d')
    ax.plot(*orbit.T,color=INK,lw=1.65)
    for fraction in (.08,.39,.71):
        j=int(fraction*len(orbit));v=orbit[j+7]-orbit[j]
        ax.quiver(*orbit[j],*v,color=COLORS[0],arrow_length_ratio=.55,linewidth=1.5)
    ax.set(xlabel=r'$X_1$',ylabel=r'$X_2$',zlabel=r'$X_5$') # One-based manuscript notation.
    ax.set_xticks([-2,-1,0,1]);ax.set_yticks([-2,-1,0,1]);ax.set_zticks([-2,-1,0,1])
    ax.set_proj_type('ortho');ax.set_box_aspect(np.ptp(orbit,axis=0));ax.view_init(24,-48)
    ax.tick_params(labelsize=11,pad=0)
    for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
        axis.set_pane_color((.965,.978,.985,1))
        axis._axinfo['grid'].update(color=(.72,.78,.83,.35),linewidth=.5)
    panel(fig,.55,.946,'B','Five cluster states')
    bx=fig.add_axes([.60,.28,.365,.54])
    for j,(row,m) in enumerate(zip([0,1,4,5,7],[2,2,1,2,1])):
        x=cycle[:,16+row] if row<7 else -np.sum(cycle[:,16:],axis=1)
        bx.plot(time,x,color=COLORS[j],lw=1.5,label=f'C{j+1} ({m} '+('rows)' if m>1 else 'row)'))
    bx.set(xlabel=r'Time / $T$',ylabel='Row displacement $X$',xlim=(0,1));clean(bx)
    handles,labels=bx.get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.78,.02),ncol=3,frameon=False,fontsize=11,
               handlelength=1.6,columnspacing=.8)
    fig.text(.275,.06,'Projection of the full 23-dimensional orbit',ha='center',color=MUTED,fontsize=11)
    save(fig,out,'figure-1-cycle',outputs)
    if args.animate:
        caption=fig.text(.27,.012,'Camera rotation only; the plotted orbit is fixed',ha='center',fontsize=8,color=MUTED)
        def rotate(frame):
            ax.view_init(24,-48+360*frame/96)
            return []
        path=out/'cycle-rotation.gif'
        FuncAnimation(fig,rotate,frames=96,interval=85).save(path,writer=PillowWriter(fps=12),dpi=110)
        outputs.append(path)
        im=Image.open(path)
        assert im.n_frames==96 and im.info.get('loop')==0
        caption.remove()
    plt.close(fig)

    fig=plt.figure(figsize=(10.8,4.9))
    panel(fig,.055,.94,'A',r'$N=8$: local attraction')
    ax=fig.add_axes([.09,.25,.36,.57])
    palette=['#007f8b','#ce6738','#7754a0','#2f8b60','#c19925','#b65179']
    for i,(h,col) in enumerate(zip(histories,palette)):
        ax.plot(h['distance_time']/period,h['full_state_orbit_distance'],color=col,lw=1.4,label=f'{i+1}')
    ax.set(yscale='log',xlabel=r'Time / $T$',ylabel='Full-state orbit distance',xlim=(0,6),ylim=(1e-8,.4))
    clean(ax)
    ax.legend(title='Declared start',ncol=6,loc='lower left',frameon=False,fontsize=11,title_fontsize=11,
              handlelength=1,columnspacing=.6)
    panel(fig,.54,.94,'B',r'$N=16$: amplified differences')
    bx=fig.add_axes([.59,.25,.365,.57])
    for i,(h,start) in enumerate(zip(new,p['starts'])):
        col=COLORS[start['cluster']]
        # Both signs and both amplitudes are plotted; signs intentionally
        # coincide in absolute value while nonlinear differences are retained.
        bx.plot(h['returns']/period,abs(h['observed'])/start['amplitude'],
                color=col,marker='o' if start['sign']>0 else 'x',ms=7 if start['cluster']==2 else 4,
                mfc='white' if start['cluster']==2 else col,
                ls='-' if start['amplitude']==1e-8 else ':',lw=1.0,alpha=.8)
    steps=np.arange(5);mu=abs(p['starts'][0]['mu'])
    bx.plot(steps,mu**steps,color=INK,lw=1.3,ls='--',label=r'Linear prediction: $|\mu|^n$')
    bx.set(yscale='log',xlabel=r'Completed return periods $n$',ylabel=r'$|a(nT)| / \varepsilon$',
           xlim=(-.1,4.1),ylim=(.55,1e6),xticks=steps)
    clean(bx)
    bx.legend(loc='lower right',frameon=False,fontsize=11)
    fig.text(.59,.065,'C3 / C5; both signs; '+r'$\varepsilon=10^{-6},10^{-8}$',fontsize=11,color=MUTED)
    fig.text(.09,.065,'Phase minimized in all 24 physical coordinates.',fontsize=10,color=MUTED)
    fig.text(.59,.025,'All eight declared trajectories are shown.',fontsize=10,color=MUTED)
    save(fig,out,'figure-2-dynamics',outputs);plt.close(fig)

    fig=plt.figure(figsize=(10.8,6.7))
    panel(fig,.045,.956,'A','The same cluster weights')
    ax=fig.add_axes([.06,.54,.46,.34]);ax.set_xlim(-.9,4.6);ax.set_ylim(-.15,1.20);ax.axis('off')
    for row,(mult,y) in enumerate([([2,2,1,2,1],.87),([4,4,2,4,2],.18)]):
        ax.text(-.8,y,rf'$N={8*(row+1)}$',ha='left',va='center',fontsize=13,fontweight='bold')
        for j,m in enumerate(mult):
            if m==1: offsets=[(0,0)]
            elif m==2: offsets=[(-.10,0),(.10,0)]
            else: offsets=[(-.10,-.08),(.10,-.08),(-.10,.08),(.10,.08)]
            for dx,dy in offsets:
                ax.scatter(j+dx,y+dy,s=85,facecolor=COLORS[j],edgecolor='white',lw=.8,zorder=3)
            if row==0:ax.text(j,1.11,f'C{j+1}',ha='center',fontsize=12,color=COLORS[j])
            if row==1 and j in (2,4):
                ax.annotate('',xy=(j+.32,y),xytext=(j-.32,y),
                            arrowprops=dict(arrowstyle='<->',color='#b63842',lw=1.5),zorder=2)
        ax.text(2,-.12 if row else .64,'4 + 4 + 2 + 4 + 2' if row else '2 + 2 + 1 + 2 + 1',
                ha='center',fontsize=12,color=MUTED)
    ax.annotate('',xy=(2,.36),xytext=(2,.57),arrowprops=dict(arrowstyle='->',color=MUTED,lw=1.4))
    panel(fig,.565,.956,'B','Certified instability count')
    bx=fig.add_axes([.63,.59,.32,.28]);ks=np.arange(1,9)
    bx.plot(8*ks,2*ks-2,color='#b63842',marker='o',lw=1.7,ms=4)
    bx.set(xlabel='Number of rows $N=8k$',ylabel='Unstable directions',xticks=[8,16,32,48,64],
           yticks=[0,2,6,10,14],ylim=(-.6,15));clean(bx)
    bx.text(.05,.82,r'$u=2k-2$',transform=bx.transAxes,fontsize=18,color='#b63842')
    fig.text(.63,.48,r'One flow direction; $22k$ stable.',fontsize=11,color=MUTED)
    panel(fig,.045,.435,'C','Sign reversal and amplification')
    cx=fig.add_axes([.09,.12,.50,.23]);steps=np.arange(5)
    for h,start in zip(new,p['starts']):
        if start['amplitude']!=1e-8:continue
        cx.plot(steps,h['observed']/start['amplitude'],color=COLORS[start['cluster']],
                marker='o' if start['sign']>0 else 'x',ms=7 if start['cluster']==2 else 4,
                mfc='white' if start['cluster']==2 else COLORS[start['cluster']],
                ls='-' if start['sign']>0 else '--',alpha=.8,lw=1.2)
    cx.axhline(0,color='#9caab4',lw=.6)
    cx.set(yscale='symlog',xlabel=r'Completed return periods $n$',ylabel=r'$a(nT) / \varepsilon$',
           xticks=steps,yticks=[-1e5,-1e2,0,1e2,1e5],ylim=(-6e5,6e5),xlim=(-.15,4.15));clean(cx)
    fig.text(.66,.32,r'$-23.084 < \mu < -23.083$',fontsize=18,color=INK)
    fig.text(.66,.245,'Certified for both new blocks',fontsize=11,color=MUTED)
    fig.text(.66,.18,'The other two multipliers satisfy',fontsize=11,color=MUTED)
    fig.text(.66,.125,r'$|\mu_{\mathrm{stable}}|<0.13$',fontsize=18,color=COLORS[0])
    save(fig,out,'figure-3-replication',outputs);plt.close(fig)
    inputs=oldpaths+newpaths+[attraction/'reference.npz',attraction/'analysis.json',nonlinear/'analysis.json',
                             nonlinear/'protocol.json',ROOT/'results/motor-replication-article/exact-index-recheck.json']
    record=dict(renderer_sha256=sha(__file__),
                input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in inputs},
                outputs={p.relative_to(ROOT).as_posix():sha(p) for p in outputs},
                qualification='All curves are numerical renderings; exact bounds and counts are supplied by separate interval certificates. Vector PDF and SVG figures retain text and paths.')
    (out/'manifest.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('Rendered',len(outputs),'assets',flush=True)


if __name__=='__main__':
    main()

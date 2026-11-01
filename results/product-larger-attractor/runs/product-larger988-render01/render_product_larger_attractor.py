"""Evidence and limitations for a larger attracting-set interpretation of Case 988."""
import argparse,hashlib,json,subprocess,sys
from fractions import Fraction
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1]
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    read=lambda p:json.loads(p.read_text())
    source=np.load(args.source);reference=source['states'];meta=read(args.metadata)
    center=np.array(meta['center']);camera=np.array(meta['basis'])@np.array(meta['projection']).T
    project=lambda x:(x-center)@camera
    probe=read(args.probe/'report.json');fan=np.load(args.probe/'unstable-fan.npz')
    long=np.load(args.probe/'branch-long.npz')['states'];rec=read(args.recurrence/'report.json')
    recurrence=np.load(args.recurrence/'recurrence.npz');trail=read(args.continuation/'trajectory.json')
    split=probe['samples_per_branch'];early=fan['states'][:,:129]
    xy=project(reference);first=project(early);late=project(long)
    # Every state in A--C fits in the common window, including rare excursions.
    half=1.04*max(np.max(abs(v)) for v in [xy,first,late])
    colors=['#a64f35','#168492'];background='#bccbd0';ink='#173c50'
    fig,axes=plt.subplots(2,3,figsize=(16,10.5));fig.subplots_adjust(left=.065,right=.985,top=.87,bottom=.14,wspace=.27,hspace=.35)
    for ax in axes[0]:
        ax.set(xlim=(-half,half),ylim=(-half,half),xlabel='View coordinate 1',ylabel='View coordinate 2')
        ax.set_aspect('equal');ax.grid(alpha=.08)
    ax=axes[0,0]
    ax.scatter(xy[::4,0],xy[::4,1],s=.18,color='#497d8a',alpha=.20,edgecolors='none',rasterized=True)
    p=project(recurrence['targets'][:4]);ax.scatter(p[:,0],p[:,1],marker='D',s=25,facecolors='white',edgecolors='#8d3b28',linewidths=1,zorder=4)
    ax.set_title('A   Reference orbit and saddle cycle',loc='left',fontsize=11.5)
    ax.text(.03,.04,'Diamonds: the four phases of cycle 14',transform=ax.transAxes,fontsize=8,color=ink)
    ax=axes[0,1];ax.scatter(xy[::10,0],xy[::10,1],s=.2,color=background,alpha=.16,edgecolors='none',rasterized=True)
    for k,color in enumerate(colors):
        z=first[k*split:(k+1)*split].reshape(-1,2)
        ax.scatter(z[:,0],z[:,1],s=.35,color=color,alpha=.45,edgecolors='none',rasterized=True)
    ax.set_title('B   Departures spread across the wings',loc='left',fontsize=11.5)
    ax.text(.03,.04,'128 launches per side · steps 88–216',transform=ax.transAxes,fontsize=8,color=ink)
    ax=axes[0,2]
    # Alternate drawing small batches to avoid making the second branch obscure
    # the first one systematically. The full million states per branch are shown.
    for offset in range(0,1000000,20000):
        for k,color in enumerate(colors):
            z=late[k,offset:offset+20000];ax.scatter(z[:,0],z[:,1],s=.12,color=color,alpha=.08,edgecolors='none',rasterized=True)
    ax.set_title('C   Both continuations recover the shape',loc='left',fontsize=11.5)
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=c,label=l) for c,l in zip(colors,['Positive launch coordinate','Negative launch coordinate'])],frameon=False,fontsize=7.8,loc='lower left')
    ax=axes[1,0]
    def decode(v):
        return [float(Fraction((-1 if s else 1)*m)*Fraction(2)**e) for s,m,e,_ in v['binary']]
    selected=[r for r in trail if r['iteration']>=80];n=np.array([r['iteration'] for r in selected]);d=np.array([decode(r['distance']) for r in selected])
    ax.fill_between(n,d[:,0],d[:,1],color=colors[0],alpha=.35)
    ax.plot(n,d[:,1],color=colors[0],marker='.',ms=3,lw=.8)
    ax.scatter([200],[d[-1,1]],color=colors[0],s=24,zorder=4)
    ax.annotate('200 steps: < 0.000011199',(200,d[-1,1]),xytext=(83,6e-6),fontsize=8,arrowprops=dict(arrowstyle='-',color='#777777'))
    ax.set(yscale='log',xlabel='Original map iterations',ylabel='Distance to a saved reference state',title='D   An exact unstable point approaches the cloud',ylim=(3e-6,.03));ax.title.set_fontsize(10.6);ax.title.set_position((.5,1.04));ax.grid(alpha=.12)
    ax=axes[1,1]
    for k,color in enumerate(colors):
        h=recurrence[f'history_{k}'][:,0];steps=np.arange(1,len(h)+1)*rec['block']
        ax.step(steps/1e6,h,where='post',color=color,lw=1.3,label=f'Independent start {k+1}')
    ax.axhline(1e-5,color='#555555',ls='--',lw=.8)
    ax.text(.98,.04,'Dashed: illustrative 0.00001 scale',transform=ax.transAxes,ha='right',fontsize=7.6,color='#555555')
    ax.set(yscale='log',xlabel='Observed iterations (millions)',ylabel='Closest distance to the saddle phase p',title='E   Reverse recurrence remains unresolved',ylim=(5e-6,.004));ax.title.set_fontsize(11);ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.12)
    ax=axes[1,2];l=probe['long_branch'];vals=[100*v['total_variation'] for v in l['grid_comparison']]+[100*l['grid_control']['total_variation']]
    ax.bar(range(3),vals,color=colors+['#83939b'],width=.6)
    for i,v in enumerate(vals):ax.text(i,v+.09,f'{v:.2f}%',ha='center',fontsize=10)
    ax.set(xticks=range(3),xticklabels=['Positive\nbranch','Negative\nbranch','Reference\nhalves'],ylim=(0,3.65),ylabel='Coarse-grid total variation (%)',title='F   Long-run occupancy matches the control');ax.title.set_fontsize(10.8);ax.grid(axis='y',alpha=.12)
    fig.suptitle('Does the unstable core belong to a larger attractor?',x=.065,ha='left',fontsize=21,color=ink)
    fig.text(.065,.925,'Case 988 · outgoing branches reproduce the wings; membership of the core remains an open question.',fontsize=11.5,color='#555555')
    fig.text(.065,.066,'A–C: native states, one fixed orthogonal camera and equal spatial units. Branch traces are numerical approximations.\nD: interval-certified continuation of one true unstable point; the target is stored numerical data. E: two fresh 50-million-step runs.\nF: cells of width 0.02; one million states per sample. These comparisons establish neither a trapping region nor invariant-set membership.',fontsize=9,color='#555555',linespacing=1.5)
    asset=root/'docs/assets/diagnostics/product-alphabet-larger-attractor.png';fig.savefig(asset,dpi=180,facecolor='white');plt.close(fig)
    subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',asset.relative_to(root).as_posix(),'--role','diagnostic'],check=True,capture_output=True)
    (args.out/asset.name).write_bytes(asset.read_bytes());(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    report=dict(status='rendered',asset=asset.relative_to(root).as_posix(),sha256=sha(asset),camera=camera.tolist(),center=center.tolist(),half_width=float(half),
        sources={str(p):sha(p) for p in [args.source,args.metadata,args.probe/'report.json',args.probe/'unstable-fan.npz',args.probe/'branch-long.npz',args.recurrence/'report.json',args.recurrence/'recurrence.npz',args.continuation/'trajectory.json']},
        scope='Diagnostic comparison, not an additional portrait or proof of attraction.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','metadata','probe','recurrence','continuation','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())

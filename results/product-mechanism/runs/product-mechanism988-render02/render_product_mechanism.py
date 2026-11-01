"""Diagnostic figures from saved states, critical images and parameter slices."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    with np.load(args.source) as f:states=f['states'];beta,tau=f['params'];levels=f['levels'][:,:2]
    with np.load(args.diagnostics) as f:values=f['values']
    with np.load(args.critical) as f:images=f['snapshots'];steps=f['checkpoints'];tags=f['tags'];distance=f['reference_distances']
    meta=json.loads(args.metadata.read_text(encoding='utf-8'));camera=np.array(meta['basis'])@np.array(meta['projection']).T;center=np.array(meta['center'])
    xy=(states-center)@camera;projected=(images-center)@camera;near=distance[0]<=.005
    half=float(max(abs(xy).max(),abs(projected[[0,1,3,5]][:,near]).max())*1.06)
    palette=np.array(['#287d8e','#bb653f','#75619b']);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(15,9.8));fig.subplots_adjust(left=.055,right=.975,top=.90,bottom=.12,hspace=.30,wspace=.26)
    for ax,n,letter in zip(axes.flat[:4],[0,1,4,32],'ABCD'):
        idx=int(np.flatnonzero(steps==n)[0]);ax.scatter(xy[::8,0],xy[::8,1],s=.08,c='#164e63',alpha=.12,edgecolors='none',rasterized=True)
        z=projected[idx,near];ax.scatter(z[:,0],z[:,1],s=1.7,c=palette[tags[near]],alpha=.65,edgecolors='none',rasterized=True)
        ax.set(xlim=(-half,half),ylim=(-half,half),xlabel='View coordinate 1',ylabel='View coordinate 2',title=f'{letter}   Critical sources' if n==0 else f'{letter}   After {n} iteration'+('' if n==1 else 's'));ax.set_aspect('equal')
    handles=[Line2D([],[],marker='o',ls='',color=c,label=f'{j}-directed root scan') for j,c in zip('xyz',palette)]
    axes[0,0].legend(handles=handles,loc='lower left',frameon=False,fontsize=8)
    ax=axes[1,1];active=(values[:,3:6]>1).sum(1);selected=np.arange(0,len(states),8)
    for mask,color,label in [(active[selected]==0,'#70a9b7','All posterior slopes at most 1'),(active[selected]>0,'#a64f35','At least one posterior slope above 1')]:
        z=xy[selected[mask]];ax.scatter(z[:,0],z[:,1],s=.16,c=color,alpha=.35,edgecolors='none',rasterized=True,label=label)
    ax.set(xlim=(-half,half),ylim=(-half,half),xlabel='View coordinate 1',ylabel='View coordinate 2',title='E   Where posterior switching occurs');ax.set_aspect('equal')
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=c,label=t) for c,t in [('#70a9b7','All posterior slopes at most 1'),('#a64f35','At least one posterior slope above 1')]],frameon=False,fontsize=7.4,loc='lower left')
    ax=axes[1,2];counts=[];means=[]
    for k in range(4):
        mask=active[1000:]==k;counts.append(int(mask.sum()));means.append(float(values[1000:,2][mask].mean()))
    ax.bar(np.arange(4),means,color=['#70a9b7','#bb653f','#a64f35','#75432e'],width=.6);ax.axhline(0,c='#777777',lw=.7)
    for k,(m,c) in enumerate(zip(means,counts)):ax.text(k,m+(.08 if m>=0 else -.09),f'{m:+.3f}\n{100*c/sum(counts):.2f}% of states',ha='center',va='bottom' if m>=0 else 'top',fontsize=8)
    ax.set(xticks=range(4),ylim=(-.65,2.35),xlabel='Number of posterior slopes above 1',ylabel='Mean one-step log tangent growth',title='F   Stretching is concentrated at switches');ax.grid(axis='y',alpha=.12)
    fig.suptitle('How the folded wings form · Case 988',x=.055,ha='left',fontsize=20,color='#122f42')
    fig.text(.055,.945,'Actual discrete states and numerical critical-set images; fixed orthogonal camera and equal spatial units.',fontsize=11,color='#555555')
    fig.text(.055,.033,'A–D: 3,310 critical-source proposals within 0.005 of the saved cloud; gray-blue background: one reference orbit.\nFinite meshes and cloud distances do not certify a boundary or attraction. E–F: two million states; first 1,000 omitted from tangent statistics.',fontsize=9,color='#555555')
    paths=[]
    def save(fig,name):
        target=root/'docs/assets/diagnostics'/name;fig.savefig(target,dpi=170,facecolor='white');plt.close(fig)
        subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',target.relative_to(root).as_posix(),'--role','diagnostic'],check=True,capture_output=True)
        (args.out/name).write_bytes(target.read_bytes());paths.append(dict(path=target.relative_to(root).as_posix(),sha256=sha(target)))
    save(fig,'product-alphabet-mechanism.png')
    with np.load(args.scan) as f:params=f['params'];tails=f['tails'];rates=f['rates'];active_counts=f['active'];periods=f['periods']
    split=len(params)//2;colors=['#164e63','#bb653f','#75619b'];fig,axes=plt.subplots(3,2,figsize=(15,10.5),sharex='col',gridspec_kw={'height_ratios':[1.4,1,.75]})
    fig.subplots_adjust(left=.065,right=.98,top=.88,bottom=.12,hspace=.15,wspace=.20)
    for col,idx in enumerate([np.arange(split),np.arange(split,len(params))]):
        parameter=params[idx,col];principal=[beta,tau][col]
        for s in range(3):
            ax=axes[0,col];ys=tails[idx,s,-512:,0];ax.scatter(np.repeat(parameter,ys.shape[1]),ys.ravel(),s=.25,c=colors[s],alpha=.45,edgecolors='none',rasterized=True)
            axes[1,col].plot(parameter,rates[idx,s].mean(1),c=colors[s],lw=.65,alpha=.9,label=f'Start {s+1}')
            fraction=1-active_counts[idx,s,0]/active_counts[idx,s].sum(1);axes[2,col].plot(parameter,fraction,c=colors[s],lw=.7)
        for row in range(3):axes[row,col].axvline(principal,c='#222222',ls='--',lw=.9);axes[row,col].grid(alpha=.12)
        if col==1:
            for row in range(3):axes[row,col].set_xscale('log')
        axes[0,col].set_title('Varying reflection weight β' if col==0 else 'Varying posterior noise τ',loc='left',fontsize=13)
        axes[2,col].set_xlabel('Reflection weight β' if col==0 else 'Posterior noise τ');axes[1,col].axhline(0,c='#777777',lw=.7)
    axes[0,0].set_ylabel('Last 512 x-coordinates');axes[1,0].set_ylabel('Finite leading growth');axes[2,0].set_ylabel('Fraction with a slope above 1')
    axes[1,1].legend(frameon=False,fontsize=9,loc='lower right')
    fig.suptitle('Geometry changes across chaotic and periodic windows',x=.065,ha='left',fontsize=20,color='#122f42')
    fig.text(.065,.925,'Case 988 alphabet fixed · 242 values per slice · three independent starts per value · dashed lines: the principal parameters',fontsize=10.5,color='#555555')
    fig.text(.065,.035,'Each cohort discards 20,000 steps and observes 80,000. Both parameters are varied separately; the other stays fixed.\nConnected diagnostic curves guide the eye between tested values. Finite-time growth and apparent periods are not bifurcation certificates.',fontsize=9,color='#555555')
    save(fig,'product-alphabet-parameter-mechanism.png')
    report=dict(status='rendered',sources={str(p):sha(p) for p in [args.source,args.diagnostics,args.critical,args.metadata,args.scan]},figures=paths,
                critical_source_selection={'distance_at_most':.005,'count':int(near.sum()),'reference_states':len(states)},camera=camera.tolist(),center=center.tolist(),half_width=half,
                conditional_means=means,conditional_counts=counts,scope=__doc__)
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());(args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(paths))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','diagnostics','critical','metadata','scan','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())

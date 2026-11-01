"""Diagnostic: finite chaotic survival and a certified singular preimage nearby."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from audit_product_orientation import interval


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1]
    read=lambda p:json.loads(p.read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    with np.load(args.source) as f:reference=f['states']
    with np.load(args.fate/'fate.npz') as f:rates=f['rates']
    meta=read(args.metadata);report=read(args.chain/'report.json');chain=read(args.chain/'chain.json');pole=read(args.singular)
    center=np.array(meta['center']);camera=np.array(meta['basis'])@np.array(meta['projection']).T
    project=lambda x:(x-center)@camera
    trail=np.vstack([np.array([[float(v) for v in row['center']] for row in chain])[::-1],np.array(pole['center'],float)])
    w=np.array(report['reference_state']);z=trail[0];xy=project(reference);path=project(trail)
    ink='#173c50';teal='#168492';orange='#a64f35';gray='#b8c8ce'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelcolor':ink,'text.color':ink,'axes.edgecolor':'#80949d','axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(14,12));fig.subplots_adjust(left=.09,right=.97,bottom=.18,top=.865,hspace=.42,wspace=.30)
    ax=axs[0,0];ax.scatter(xy[::4,0],xy[::4,1],s=.12,c=gray,alpha=.22,edgecolors='none',rasterized=True)
    ax.scatter(path[:-1,0],path[:-1,1],s=22,color=orange,edgecolors='white',linewidths=.4,zorder=3)
    ax.scatter(path[-1,0],path[-1,1],s=105,marker='X',color=ink,edgecolors='white',linewidths=.7,zorder=4)
    ax.scatter(path[0,0],path[0,1],s=82,facecolors='none',edgecolors=orange,lw=1.4,zorder=5)
    ax.annotate('z: 27 steps before the pole',path[0],xytext=(.02,.97),textcoords='axes fraction',va='top',fontsize=9,arrowprops=dict(arrowstyle='-',color=orange))
    ax.annotate('s: map undefined',path[-1],xytext=(.49,.06),textcoords='axes fraction',fontsize=9,arrowprops=dict(arrowstyle='-',color=ink))
    half=1.05*max(abs(xy).max(),abs(path).max());ax.set(xlim=(-half,half),ylim=(-half,half),aspect='equal',xlabel='View coordinate 1',ylabel='View coordinate 2')
    ax.set_title('A   A certified chain beside the numerical wings',loc='left',fontsize=11)
    ax=axs[0,1];delta=(z-w)/1e-5
    ax.add_patch(Circle((0,0),6,facecolor=teal,alpha=.05,edgecolor='none'));ax.add_patch(Circle((0,0),6,fill=False,edgecolor=teal,ls='--',lw=1.2))
    ax.scatter([0],[0],s=70,color=teal,zorder=3);ax.scatter([delta[0]],[delta[1]],s=80,marker='X',color=orange,zorder=3)
    ax.annotate('w: saved numerical state',(0,0),xytext=(-5.7,-2),fontsize=9,arrowprops=dict(arrowstyle='-',color=teal))
    ax.annotate('z: singular preimage',delta[:2],xytext=(-5.7,3.3),fontsize=9,arrowprops=dict(arrowstyle='-',color=orange))
    ax.text(.04,.04,'Full 3D distance < 0.000060\nThird-coordinate difference: 0.00000400',transform=ax.transAxes,fontsize=9,linespacing=1.5)
    ax.set(xlim=(-7,7),ylim=(-7,7),aspect='equal',xlabel='(x − wₓ) / 0.00001',ylabel='(y − wᵧ) / 0.00001')
    ax.set_title('B   This entire ball cannot lie in a smooth trap',loc='left',fontsize=11)
    ax=axs[1,0]
    radius=np.array([float(interval(row['normalization'])[0]) for row in chain])[::-1]
    ax.plot(np.arange(27),radius,'o-',color=orange,lw=1.1,ms=3)
    ax.scatter([27],[0],s=85,marker='X',color=ink,zorder=3);ax.axhline(0,color=gray,lw=.7)
    ax.text(.04,.05,'27 defined steps reach s.\nThe next step does not exist.',transform=ax.transAxes,fontsize=10,linespacing=1.5)
    ax.set(xlabel='Iteration from z',ylabel='Normalization radius',xlim=(-.7,28),ylim=(-.18,max(radius)*1.15),xticks=[0,5,10,15,20,25,27])
    ax.set_title('C   Exact obstruction, verified by two expressions',loc='left',fontsize=11)
    ax=axs[1,1];labels=['Global cube','Near cloud\n0.00001','Near cloud\n0.01','Near cycle\n0.01']
    for g in range(4):
        group=rates[64*g:64*(g+1)];means=group.mean(1)
        ax.vlines(g,group.min(),group.max(),color=gray,lw=9,alpha=.7)
        ax.vlines(g,means.min(),means.max(),color=teal,lw=6)
        ax.scatter(g+np.linspace(-.20,.20,64),means,s=9,color=teal,alpha=.65,edgecolors='none')
    ax.text(.04,.92,'256 starts × 4 million observed steps\nAll 5,120 growth blocks are positive',transform=ax.transAxes,va='top',fontsize=10,linespacing=1.5)
    ax.set(xlim=(-.55,3.55),xticks=range(4),xticklabels=labels,ylabel='Leading finite-time growth per step',ylim=(.238,.268))
    ax.set_title('D   Persistent finite chaos; no attraction theorem',loc='left',fontsize=11)
    ax.text(.04,.04,'Teal: trajectory means · gray: 200,000-step blocks',transform=ax.transAxes,fontsize=8.5)
    for ax in axs.flat:ax.grid(alpha=.10)
    fig.suptitle('Case 988: is the four-wing regime an attractor?',x=.09,ha='left',fontsize=21,color=ink)
    fig.text(.09,.922,'Strong finite evidence of chaos. A certified singular preimage constrains any proposed trapping neighborhood.',fontsize=11,color='#526a76')
    fig.text(.09,.055,'A: fixed orthogonal view; discrete states are not connected by a flow curve. B: native x–y projection of a 3D ball of radius 0.000060.\nC: interval lower bounds (widths are below plotting resolution); the final singularity is certified separately.\nD: 20,000 warm-up steps per start; no tested periodic capture. This finite search does not rule out long transients.\nThe reference point w is stored numerical data. The obstruction does not exclude every smaller trap or a Milnor attractor.',fontsize=9,color='#526a76',linespacing=1.55)
    asset=root/'docs/assets/diagnostics/product-alphabet-attractor-status.png';fig.savefig(asset,dpi=180,facecolor='white');plt.close(fig)
    subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',asset.relative_to(root).as_posix(),'--role','diagnostic'],check=True,capture_output=True)
    (args.out/asset.name).write_bytes(asset.read_bytes());(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    result=dict(status='rendered',asset=asset.relative_to(root).as_posix(),sha256=sha(asset),camera=camera.tolist(),center=center.tolist(),half_width=float(half),
        native_zoom_difference=delta.tolist(),normalization_lower_bounds=radius.tolist(),
        sources={str(p):sha(p) for p in [args.source,args.metadata,args.fate/'fate.npz',args.chain/'chain.json',args.chain/'report.json',args.singular]},scope=__doc__)
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(dict(status=result['status'],sha256=result['sha256'])))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','metadata','fate','chain','singular','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())

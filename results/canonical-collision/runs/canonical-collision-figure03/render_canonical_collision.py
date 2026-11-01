"""Diagnostic of validated discrete histories and their common stable return."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def interval(text):
    return np.array([float(v.strip()) for v in text.strip('[]').split(',')])


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    args.image.parent.mkdir(parents=True,exist_ok=True)
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    c=json.loads(args.collision.read_text(encoding='utf-8'))
    h=json.loads(args.homoclinic.read_text(encoding='utf-8'))
    s,t=map(float,c['fixed_point']['center']);p=np.array([s,s,t])
    excursion=np.array([[float(v) for v in row['center']] for row in h['crossing']['last40']])
    reflected=excursion[:,[1,0,2]]
    crossing=np.array([interval(v).mean() for v in h['crossing']['box']])
    boxes=[h['crossing']['box']]+[row['box'] for row in h['return_to_saddle']['trail']]
    returning=np.array([[interval(v).mean() for v in row] for row in boxes])
    basis=np.array([[1,-1,0],[1,1,-2]],float)/np.array([[np.sqrt(2)],[np.sqrt(6)]])
    project=lambda x:x@basis.T
    a,b,g,P,C=map(project,[excursion,reflected,returning,p,crossing])
    colors=dict(ink='#233f4b',muted='#71848b',orange='#c66c31',blue='#2b7c9a',green='#5a8767')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelcolor':colors['ink'],'text.color':colors['ink'],'axes.edgecolor':'#cbd4d7','xtick.color':colors['muted'],'ytick.color':colors['muted']})
    fig=plt.figure(figsize=(14,6.2),facecolor='white')
    grid=fig.add_gridspec(1,2,width_ratios=[1.45,1],left=.075,right=.97,top=.82,bottom=.19,wspace=.26)
    ax=fig.add_subplot(grid[0]);bx=fig.add_subplot(grid[1])
    for q in [ax,bx]:
        q.spines[['top','right']].set_visible(False)
        q.grid(color='#e9eef0',linewidth=.65,alpha=.8)
        q.set_axisbelow(True)
    ax.axvline(0,color='#cbd4d7',linewidth=1,zorder=1)
    for data,color,label in [(a,colors['orange'],'One past'),(b,colors['blue'],'Reflected past')]:
        ax.plot(data[:,0],data[:,1],color=color,linewidth=.8,linestyle=(0,(3,3)),alpha=.65,zorder=2)
        ax.scatter(data[:,0],data[:,1],s=20,c=color,edgecolors='white',linewidths=.35,label=label,zorder=3)
    ax.plot(g[:,0],g[:,1],color=colors['green'],linewidth=2,zorder=4)
    ax.scatter(g[::3,0],g[::3,1],color=colors['green'],s=14,label='Common future',zorder=5)
    ax.scatter(*P,s=180,marker='*',color=colors['ink'],edgecolors='white',linewidths=.7,zorder=6)
    ax.scatter(*C,s=65,marker='D',color=colors['green'],edgecolors='white',linewidths=.8,zorder=7)
    ax.annotate('Saddle p',xy=P,xytext=(8,9),textcoords='offset points',fontsize=10,fontweight='bold')
    ax.set_xlabel(r'$(x-y)/\sqrt{2}$');ax.set_ylabel(r'$(x+y-2z)/\sqrt{6}$')
    ax.set_aspect('equal',adjustable='box');ax.set_anchor('N');ax.margins(.12)
    ax.set_title('A   Two histories merge by symmetry',loc='left',fontsize=13,pad=17)
    ax.legend(loc='upper left',frameon=False,fontsize=9,handletextpad=.5)
    inverse=np.array([[float(Fraction(v)) for v in row] for row in c['unstable_graph']['inverse_fraction']])
    distances=np.max(np.abs((returning-p)@inverse.T)[:,1:],axis=1)
    bx.semilogy(np.arange(len(distances)),distances,color=colors['green'],linewidth=2)
    bx.scatter(np.arange(0,len(distances),5),distances[::5],s=20,color=colors['green'],edgecolors='white',linewidths=.5)
    radius=float(h['return_to_saddle']['local_radius']);n=h['return_to_saddle']['iterations']
    bx.axhline(radius,color=colors['muted'],linestyle='--',linewidth=1)
    bx.annotate('Certified contracting square',xy=(n,radius),xytext=(.04,.14),textcoords='axes fraction',fontsize=10,color=colors['ink'])
    bx.set_xlabel('Steps after the merger');bx.set_ylabel('Distance in stable coordinates')
    bx.set_xlim(0,n+3);bx.set_ylim(radius/2,max(distances)*1.6)
    bx.set_title('B   The common future returns to p',loc='left',fontsize=13,pad=17)
    fig.text(.075,.94,'A certified homoclinic connection',fontsize=22,fontweight='normal',color=colors['ink'])
    fig.text(.075,.89,'Canonical cube · Case 653 · discrete posterior-reflection map',fontsize=12,color=colors['muted'])
    fig.text(.075,.30,'421 steps to the merger · 61 more to local contraction',fontsize=12,color=colors['green'])
    fig.text(.075,.25,'The diamond marks the shared state; the star marks the saddle.',fontsize=10,color=colors['muted'])
    fig.text(.075,.09,'Markers show enclosure centers: the last 40 outgoing steps and the return. Dashed segments are discrete-iteration guides.',fontsize=10,color=colors['muted'])
    fig.text(.075,.055,'The two pasts are distinct. Their merged future enters a contracting neighborhood of the same saddle.',fontsize=10,color=colors['muted'])
    fig.savefig(args.image,dpi=180,facecolor='white');plt.close(fig)
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    report=dict(status='rendered',collision_sha256=sha(args.collision),homoclinic_sha256=sha(args.homoclinic),image=str(args.image),image_sha256=sha(args.image),projection=basis.tolist(),scope='Diagnostic of the certified homoclinic enclosure centers, not a portrait of the full attracting set. Two-dimensional orthogonal projection with equal units; dashed lines are iteration guides. The certificate, not pixel placement or rounded plotting arithmetic, establishes the connection.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['collision','homoclinic','out','image']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())

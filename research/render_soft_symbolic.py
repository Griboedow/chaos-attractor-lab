"""Actual-state cycle geometry and a schematic certified covering graph."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from render_portraits import COLOR
INK='#203d4b';GOLD='#c17b23';PURPLE='#8358a1';BLUE='#167b9a'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);args.figures.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':'#b8cbd2','xtick.color':INK,'ytick.color':INK})
    with np.load(args.states) as f:states=f['states'];params=f['params'];points=f['points']
    view=json.loads(args.view.read_text());center=np.array(view['center']);rotation=np.array(view['basis'])@np.array(view['projection']).T;assert np.allclose(rotation.T@rotation,np.eye(2),rtol=0,atol=1e-14);xy=(states-center)@rotation
    cycles=json.loads((args.cycles/'report.json').read_text());assert params.tolist()==cycles['params'] and points.tolist()==cycles['points'];orbits=[]
    for row in cycles['cases']:
        orbit=np.array([[sum(map(float,v.strip('[]').split(',')))/2 for v in q['box']] for q in row['trapping_membership']]);orbits.append(orbit)
    allxy=np.vstack([xy]+[(o-center)@rotation for o in orbits]);lo=allxy.min(0);hi=allxy.max(0);span=hi-lo;lo-=.06*span;hi+=.06*span;bins=(1500,int(np.ceil(1500*(hi[1]-lo[1])/(hi[0]-lo[0]))));hist,xe,ye=np.histogram2d(xy[:,0],xy[:,1],bins=bins,range=list(zip(lo,hi)));assert hist.sum()==len(states);level=np.log1p(hist);vmax=np.quantile(level[hist>0],.9995);intensity=np.minimum(level/vmax,1)**.72
    fig=plt.figure(figsize=(13.2,5.9));ax=fig.add_axes([.035,.22,.61,.57]);ax.imshow(intensity.T,origin='lower',extent=[lo[0],hi[0],lo[1],hi[1]],cmap=COLOR,vmin=0,vmax=1,interpolation='nearest');ax.set_aspect('equal');ax.axis('off');labels=['P: fixed point, index 1','A: period 5, index 2','B: period 5, index 2'];colors=[INK,GOLD,PURPLE];markers=['D','o','s']
    for orbit,label,color,marker in zip(orbits,labels,colors,markers):
        q=(orbit-center)@rotation;ax.scatter(q[:,0],q[:,1],s=42,c=color,marker=marker,edgecolors='white',linewidths=.8,zorder=3,label=label)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.08),ncol=1,frameon=False,fontsize=10)
    chart=fig.add_axes([.73,.27,.24,.46]);chart.axvline(1,color='#b2bdc3',lw=1.1,ls='--');chart.axvspan(1,10,color='#f5ead9',alpha=.65)
    for row,color,marker,y in zip(cycles['cases'],colors,markers,[2,1,0]):
        multipliers=abs(np.array(row['numerical_multipliers'],float));chart.scatter(multipliers,np.full(3,y),c=color,marker=marker,s=55,zorder=3,edgecolors='white',linewidths=.5)
    chart.set_xscale('log');chart.set_xlim(1e-7,12);chart.set_ylim(-.5,2.5);chart.set_yticks([2,1,0],['P','A','B']);chart.set_xticks([1e-6,1e-3,1]);chart.set_xlabel('Absolute return multiplier',labelpad=12);chart.spines[['top','right','left']].set_visible(False);chart.tick_params(axis='y',length=0);chart.grid(axis='x',alpha=.15);chart.set_title('One or two expanding directions',fontsize=11,pad=19)
    fig.text(.04,.935,'Periodic structure inside a folded band',fontsize=19,weight='bold');fig.text(.04,.88,'Case 156  /  2,000,000 consecutive states  /  certified periodic points',fontsize=11,color='#55727f');fig.text(.73,.09,'Right of 1: expansion over one\ncomplete return of that orbit.',fontsize=10,color='#55727f');fig.text(.04,.035,'Orthographic projection with equal spatial units. Markers are discrete states; no connecting orbit curves.',fontsize=10,color='#55727f')
    geometry=args.figures/'soft-periodic-structure.png';fig.savefig(geometry,dpi=200,facecolor='white');plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,5.5));fig.subplots_adjust(left=.03,right=.98,bottom=.09,top=.85);ax.set_xlim(-3.1,3.1);ax.set_ylim(-1.6,1.6);ax.set_aspect('equal');ax.axis('off');angles=np.linspace(np.pi,-np.pi,7)[:-1];pos=np.c_[1.05*np.cos(angles),1.05*np.sin(angles)]
    for i in range(6):ax.add_patch(FancyArrowPatch(pos[i],pos[(i+1)%6],connectionstyle='arc3,rad=-.08',arrowstyle='-|>',mutation_scale=16,color=BLUE,lw=2,shrinkA=16,shrinkB=16))
    x,y=pos[0];ax.add_patch(FancyArrowPatch((x-.13,y+.14),(x-.13,y-.14),connectionstyle='arc3,rad=2.8',arrowstyle='-|>',mutation_scale=16,color=GOLD,lw=2))
    for i,(x,y) in enumerate(pos):ax.add_patch(Circle((x,y),.17,facecolor=INK if i==0 else '#eef5f7',edgecolor=INK,lw=1));ax.text(x,y,f'$N_{i}$',ha='center',va='center',color='white' if i==0 else INK,fontsize=12)
    ax.text(0,.13,'Each arrow',ha='center',fontsize=11);ax.text(0,-.2,'one F iteration',ha='center',fontsize=13,weight='bold');ax.text(-2.95,1.03,'0  /  STAY',color=GOLD,weight='bold');ax.text(-2.95,.71,'Six self-loops\nat N0.',linespacing=1.5);ax.text(1.65,1.03,'1  /  EXCURSION',color=BLUE,weight='bold');ax.text(1.65,.71,'Visit all six\nwindows in order.',linespacing=1.5);ax.text(0,-1.45,'Both routes take six steps and remain in the certified attracting set.',ha='center',fontsize=11)
    fig.text(.04,.94,'Two routes, arbitrary binary sequences',fontsize=19,weight='bold');fig.text(.04,.885,'Seven covering relations verified by two interval implementations. Node positions are schematic.',fontsize=10,color='#55727f');graph=args.figures/'soft-covering-graph.png';fig.savefig(graph,dpi=200,facecolor='white');plt.close(fig)
    np.savez_compressed(args.out/'display-coordinates.npz',center=center,rotation=rotation,**{f'orbit_{r["id"]}':o for r,o in zip(cycles['cases'],orbits)})
    report=dict(status='passed',states_sha256=sha(args.states),view_sha256=sha(args.view),cycle_report_sha256=sha(args.cycles/'report.json'),states_count=len(states),histogram_count=int(hist.sum()),camera=rotation.tolist(),center=center.tolist(),spatial_clipping=0,figures={p.name:sha(p) for p in [geometry,graph]},source_sha256=sha(Path(__file__)),shared_palette_source_sha256=sha(ROOT/'tools/render_portraits.py'),scope='Actual numerical states and rounded centers of certified periodic-state boxes; no interpolation, replication or synthetic geometry. The multipliers are rounded estimates backed by interval brackets. The covering graph is schematic.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps(report))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--states',type=Path,required=True);p.add_argument('--view',type=Path,required=True);p.add_argument('--cycles',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--figures',type=Path,required=True);main(p.parse_args())

"""Render actual RAAR states, a covering graph, and numerical coded itineraries."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from scipy.spatial import cKDTree
from render_twist import camera

INK='#203d4b'; BLUE='#167b9a'; GOLD='#db8a24'; PALE='#eef5f7'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main(args):
    repo=Path(__file__).resolve().parents[1]; base=args.base; out=args.out
    out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':'#b8cbd2','xtick.color':INK,'ytick.color':INK,'pdf.fonttype':42})
    reference=repo/'results/raar-circles/render/case-1051-animation.npz'
    states=np.load(reference)['states']
    view=next(r for r in json.loads((repo/'results/raar-circles/render/analysis.json').read_text())['runs'] if r['id']==1051)
    center=np.array(view['center']); rotation=np.array(view['basis'])@camera(view['azimuth'],view['elevation']).T
    cycle49=json.loads((base/'raar-symbolic-cycle-audit01/report.json').read_text())
    orbit49=np.array(cycle49['orbit_decimal'],float)
    c7=json.loads((repo/'results/raar-cycle-structure/cycle-proof/base-unstable-period-7.json').read_text())
    # Existing certificate stores decimal point orbit for display, separate from interval proof.
    if 'orbit_decimal' in c7: orbit7=np.array(c7['orbit_decimal'],float)
    else:
        from continue_raar_cycles import iterate
        orbit7=iterate(np.array(c7['center_decimal'],float),np.array(c7['parameters']),7)[3]
    xy=(states-center)@rotation; q49=(orbit49-center)@rotation; q7=(orbit7-center)@rotation
    lo=xy[:,:2].min(0);hi=xy[:,:2].max(0);span=hi-lo;lo-=.04*span;hi+=.04*span
    hist,xe,ye=np.histogram2d(xy[:,0],xy[:,1],bins=(1800,1150),range=list(zip(lo,hi)))
    assert int(hist.sum())==len(states)
    fig,ax=plt.subplots(figsize=(11,6.6));fig.subplots_adjust(left=.025,right=.98,bottom=.055,top=.89)
    palette=LinearSegmentedColormap.from_list('arch',['#ffffff','#b9dce7','#3e99b4','#0e5273','#102f41'])
    ax.imshow(np.log1p(hist.T),origin='lower',extent=[lo[0],hi[0],lo[1],hi[1]],cmap=palette,interpolation='nearest',norm=PowerNorm(.7,vmin=0,vmax=np.log1p(hist.max())),rasterized=True)
    ax.scatter(q49[:,0],q49[:,1],s=21,c=GOLD,edgecolors='white',linewidths=.45,zorder=3,label='Certified least-period-49 saddle')
    ax.scatter(q7[:,0],q7[:,1],s=40,facecolors='none',edgecolors=INK,linewidths=.9,zorder=4,label='Certified least-period-7 saddle')
    ax.set_aspect('equal');ax.set_xlim(lo[0],hi[0]);ax.set_ylim(lo[1],hi[1]);ax.axis('off')
    ax.legend(loc='lower left',frameon=True,facecolor='white',edgecolor='none',fontsize=9)
    fig.text(.04,.95,'A SYMBOLIC SKELETON INSIDE THE SERRATED ARCH',weight='bold',fontsize=17,color=INK)
    fig.text(.04,.91,'1,000,000 reference states  /  native coordinates  /  discrete points, no connecting curves',fontsize=10,color='#55727f')
    for ext in ['png','pdf']:fig.savefig(out/('raar-symbolic-geometry.'+ext),dpi=230,facecolor='white')
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(11,5.7));fig.subplots_adjust(left=.03,right=.98,bottom=.06,top=.88);ax.set_xlim(-2.4,4.5);ax.set_ylim(-1.7,1.7);ax.set_aspect('equal');ax.axis('off')
    angles=np.linspace(np.pi, -np.pi,8)[:-1];positions=np.column_stack([1.15+1.35*np.cos(angles),1.15*np.sin(angles)])
    def arrow(a,b,rad=0):ax.add_patch(FancyArrowPatch(a,b,connectionstyle=f'arc3,rad={rad}',arrowstyle='-|>',mutation_scale=15,color=BLUE,lw=2,shrinkA=17,shrinkB=17))
    for i in range(7):arrow(positions[i],positions[(i+1)%7],-.08)
    x,y=positions[0];ax.add_patch(FancyArrowPatch((x-.13,y+.15),(x-.13,y-.15),connectionstyle='arc3,rad=2.8',arrowstyle='-|>',mutation_scale=15,color=GOLD,lw=2))
    for i,(x,y) in enumerate(positions):
        ax.add_patch(Circle((x,y),.18,facecolor=INK if i==0 else PALE,edgecolor=INK,lw=1.2));ax.text(x,y,f'$N_{i}$',ha='center',va='center',color='white' if i==0 else INK,fontsize=12)
    ax.text(1.15,.1,'Every arrow',ha='center',fontsize=11);ax.text(1.15,-.19,r'$G=F^7$',ha='center',fontsize=19,weight='bold')
    ax.text(-1.9,1.15,'0  /  STAY',weight='bold',color=GOLD);ax.text(-1.9,.82,'Take the self-loop\nseven times.',fontsize=10,linespacing=1.55)
    ax.text(2.9,1.15,'1  /  EXCURSION',weight='bold',color=BLUE);ax.text(2.9,.82,'Visit all seven\nwindows in order.',fontsize=10,linespacing=1.55)
    ax.text(1.15,-1.63,'Both choices return to N0 after 49 original iterations.',ha='center',fontsize=11)
    fig.text(.04,.95,'TWO ROUTES. ANY BINARY SEQUENCE.',weight='bold',fontsize=17)
    fig.text(.04,.90,'Eight interval-certified covering relations; arrows describe a graph, not spatial orbit segments.',fontsize=10,color='#55727f')
    for ext in ['png','pdf']:fig.savefig(out/('raar-symbolic-graph.'+ext),dpi=220,facecolor='white')
    plt.close(fig)

    words=['1','001','001011']; fig,axes=plt.subplots(3,1,figsize=(11,6.3),sharey=True);fig.subplots_adjust(left=.08,right=.97,top=.84,bottom=.10,hspace=.47)
    stored={};tree=cKDTree(orbit7)
    for ax,word in zip(axes,words):
        data=np.load(base/'raar-symbolic-orbits01'/f'word-{word}.npz');orb=data['orbit'];dist=tree.query(orb)[0];stored['word_'+word]=dist
        for k,letter in enumerate(word):
            ax.axvspan(49*k-.5,49*(k+1)-.5,color=GOLD if letter=='0' else BLUE,alpha=.085)
            ax.text(49*k+24,1.035,letter,transform=ax.get_xaxis_transform(),ha='center',va='bottom',weight='bold',color=GOLD if letter=='0' else BLUE)
        assert dist.min()>1e-4 and dist.max()<.7
        ax.scatter(np.arange(len(dist)),dist,s=8,color=INK,linewidths=0);ax.set_yscale('log');ax.set_ylim(1e-4,.7);ax.set_xlim(-.5,len(dist)-.5);ax.set_xticks(np.arange(0,len(dist)+1,49));ax.set_ylabel('Distance to\n7-cycle',fontsize=9);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.text(.995,1.06,f'{len(dist)} states',transform=ax.transAxes,ha='right',fontsize=8,color='#55727f')
    axes[-1].set_xlabel('Original iteration index (one symbol = 49 steps)')
    fig.text(.06,.95,'CODED EXCURSIONS, SEEN IN ACTUAL STATES',weight='bold',fontsize=17)
    fig.text(.06,.90,'Multiple-shooting illustrations. Distances are to all seven cycle points, not to one selected phase.',fontsize=10,color='#55727f')
    for ext in ['png','pdf']:fig.savefig(out/('raar-symbolic-itineraries.'+ext),dpi=220,facecolor='white')
    plt.close(fig)
    np.savez_compressed(out/'figure-data.npz',orbit7=orbit7,orbit49=orbit49,rotation=rotation,center=center,histogram=hist.astype(np.uint32),xedges=xe,yedges=ye,**stored)
    result=dict(reference_sha256=sha(reference),reference_points=len(states),histogram_count=int(hist.sum()),cycle49_points=len(orbit49),cycle7_points=len(orbit7),words=words,source_sha256=sha(Path(__file__)),scope='Rigid projection of actual discrete states; no interpolated trajectories. The reference is finite and cycle membership in an attracting set is not certified. Graph layout is schematic. Longer coded cycles are numerical illustrations; their ill-conditioned binary64 small multipliers are not used.')
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

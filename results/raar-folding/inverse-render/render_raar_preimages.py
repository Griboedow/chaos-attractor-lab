"""Two rigid views of four certified preimages; arrows denote one-step relations."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from render_twist import camera


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);cases=json.loads(args.proof.read_text())['cases'];case=next(r for r in cases if r['label']=='target-254999');x=np.array([[float(v) for v in r['center_decimal']] for r in case['preimages']])
    with np.load(args.reference) as d:z=d['states'];target=z[254999]
    r=json.loads(args.camera.read_text())['runs'][0];center=np.array(r['center']);basis=np.array(r['basis']);colors=['#e5a330','#a45c9d','#208a52','#986b2e'];views=[(3.,.4),(.5,.5)];projected=[]
    for az,el in views:
        transform=basis@camera(az,el).T;projected.append(((z-center)@transform,(x-center)@transform,(target-center)@transform))
    half=max(np.max(np.abs(np.r_[cloud,points,y[None,:]])) for cloud,points,y in projected)*1.04
    fig,axes=plt.subplots(1,2,figsize=(10,4.6))
    for ax,(cloud,points,y),view in zip(axes,projected,views):
        ax.scatter(cloud[:,0],cloud[:,1],s=.015,alpha=.2,c='#417b8a',linewidths=0,rasterized=True)
        for j,(point,color) in enumerate(zip(points,colors)):
            ax.annotate('',xy=y,xytext=point,arrowprops=dict(arrowstyle='->',lw=1,color=color,alpha=.8));ax.scatter(*point,s=30,color=color,zorder=5,edgecolors='white',lw=.5);ax.annotate(f'x{j+1}',point,xytext=[(5,8),(5,8),(-25,-18),(16,24)][j],textcoords='offset points',fontsize=9,color=color,arrowprops=dict(arrowstyle='-',lw=.6,color=color))
        ax.scatter(*y,s=60,c='#be3636',marker='*',zorder=6);ax.annotate('same y',y,xytext=(6,-14),textcoords='offset points',fontsize=9,color='#9b2929');ax.set(xlim=(-half,half),ylim=(-half,half),aspect='equal');ax.set_axis_off()
    axes[0].set_title('Four distinct certified inputs',fontsize=12);axes[1].set_title('Same points, second rigid view',fontsize=12)
    fig.suptitle('F(x₁) = F(x₂) = F(x₃) = F(x₄) = y',fontsize=16);fig.text(.5,.045,'Arrows show a discrete update, not continuous trajectories. Gray: finite reference cloud.',ha='center',fontsize=10,color='#51666c');fig.subplots_adjust(left=.01,right=.99,bottom=.09,top=.88,wspace=.03)
    fig.savefig(args.out/'raar-four-preimages.png',dpi=250);fig.savefig(args.out/'raar-four-preimages.pdf',dpi=300);plt.close(fig)
    np.savez_compressed(args.out/'plotted-preimages.npz',preimages=x,target=target)
    row=dict(proof=str(args.proof),proof_sha256=hashlib.sha256(args.proof.read_bytes()).hexdigest(),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),camera_sha256=hashlib.sha256(args.camera.read_bytes()).hexdigest(),target_index=254999,views=views,common_halfwidth=half,reference_points=len(z),scope='Actual coordinate centers of certified boxes, two rigid orthographic cameras, common equal-unit scale. The other preimages are away from the finite reference sample; no assertion that all lie on an attracting set.')
    (args.out/'analysis.json').write_text(json.dumps(row,indent=2),encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--proof',type=Path,required=True);p.add_argument('--reference',type=Path,default=Path('results/raar-circles/render/case-1051-animation.npz'));p.add_argument('--camera',type=Path,default=Path('results/raar-circles/render/analysis.json'));p.add_argument('--out',type=Path,required=True);main(p.parse_args())

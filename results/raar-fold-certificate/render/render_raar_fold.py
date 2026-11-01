"""Actual determinant-sign densities and a certified fold, before and after one update."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from render_twist import camera
from verify_raar_folding import sha

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);proof=json.loads(args.proof.read_text());report=json.loads(args.proposal.read_text())
    with np.load(args.source) as d:z=d['states']
    with np.load(args.proposal.parent/'determinants.npz') as d:det=d['determinants'][:-1]
    cam=next(r for r in json.loads(args.camera.read_text())['runs'] if r['id']==1051);center=np.array(cam['center']);projection=np.array(cam['basis'])@camera(cam['azimuth'],cam['elevation']).T;xy=(z-center)@projection
    point=np.array([float(v) for v in proof['state_center']]);image=np.array([float(v.strip('[]').split(',')[0]) for v in proof['image_box']]);stars=[(point-center)@projection,(image-center)@projection];xb=[xy[:,0].min()-.03,xy[:,0].max()+.03];yb=[xy[:,1].min()-.03,xy[:,1].max()+.03];hist=[];blue=np.array([.18,.62,.85]);gold=np.array([1.,.64,.22]);fig,axes=plt.subplots(1,2,figsize=(11,4.8),facecolor='#071320')
    for ax,q,star,title in zip(axes,[xy[:-1],xy[1:]],stars,['Input states x','Their next states F(x)']):
        pair=[]
        for sign in [-1,1]:
            points=q[det*sign>0];h,_,_=np.histogram2d(points[:,0],points[:,1],bins=(1000,700),range=[xb,yb]);pair.append(h.T[::-1].astype(np.uint32))
        pair=np.array(pair);hist.append(pair);total=pair.sum(0);assert total.sum()==999999;w=np.divide(pair[1],total,out=np.zeros(total.shape),where=total!=0);rgb=((1-w[:,:,None])*blue+w[:,:,None]*gold)*np.clip(np.log1p(total)/np.log1p(50),0,1)[:,:,None]**.6
        ax.imshow(rgb,extent=[*xb,*yb],aspect='equal');ax.scatter(*star,c='white',marker='*',s=50,edgecolors='#071320',lw=.4,zorder=5);ax.annotate('certified fold' if title.startswith('Input') else 'its image',star,xytext=(12,-22),textcoords='offset points',color='white',fontsize=9,arrowprops=dict(arrowstyle='-',color='white',lw=.6));ax.set_title(title,color='white',fontsize=12);ax.set_axis_off()
    fig.suptitle('A genuine fold of the full three-dimensional map',color='white',fontsize=16);fig.text(.5,.065,'Blue: det DF(x) < 0     Gold: det DF(x) > 0     Mixed color: overlapping projected densities',ha='center',color='#d3e3ed',fontsize=10);fig.text(.5,.025,'Same rigid camera and spatial scale; 999,999 actual consecutive-state pairs. Fold membership in the attractor is unproved.',ha='center',color='#a7bbcc',fontsize=8);fig.subplots_adjust(left=.015,right=.985,bottom=.14,top=.88,wspace=.04)
    fig.savefig(args.out/'raar-certified-fold.png',dpi=250,facecolor=fig.get_facecolor());fig.savefig(args.out/'raar-certified-fold.pdf',dpi=300,facecolor=fig.get_facecolor());plt.close(fig);np.savez_compressed(args.out/'sign-histograms.npz',histograms=hist,projection=projection,center=center,fold_state=point,fold_image=image,x_bounds=xb,y_bounds=yb)
    meta=dict(proof_sha256=sha(args.proof),proposal_sha256=sha(args.proposal),source_sha256=sha(args.source),point_pairs=999999,scope='Density of actual input/output pairs, colored by the determinant sign at input. Star coordinates are rounded certificate centers. Mixed colors indicate a projection overlap, not proof of equal full-state preimages. No causal or attractor-membership claim.')
    (args.out/'analysis.json').write_text(json.dumps(meta,indent=2),encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--proof',type=Path,required=True);p.add_argument('--proposal',type=Path,required=True);p.add_argument('--source',type=Path,default=Path('results/raar-circles/render/case-1051-animation.npz'));p.add_argument('--camera',type=Path,default=Path('results/raar-circles/render/analysis.json'));p.add_argument('--out',type=Path,required=True);main(p.parse_args())

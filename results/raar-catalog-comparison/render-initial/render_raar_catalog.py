"""Show sampled geometric matches without treating them as equivalence tests."""
import argparse,json,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from compare_attractor_shapes import normalize
from compare_affine_shapes import whiten
from render_twist import camera

def main(args):
    repo=Path(__file__).resolve().parents[1];out=args.out;out.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.source/'report.json').read_text());saved=np.load(args.source/'aligned-samples.npz')
    with np.load(repo/'results/raar-circles/render/case-1051-animation.npz') as f:z=f['states']
    rng=np.random.default_rng(2026110122);plots={};info=[];projection=camera(.6,.35).T
    selection=[('similarity',0),('similarity',1),('affine',0)]
    for col,(kind,rank) in enumerate(selection):
        transform=normalize if kind=='similarity' else whiten
        row=report['representations'][kind]['nearest'][rank];name=row['reference'];key=f'{kind}--{name}'
        a=transform(z);b=transform(np.load(repo/f'results/attractor-novelty/reference/trajectories/{name}.npz')['points'])
        ai=saved[key+'--candidate_indices'];bi=saved[key+'--reference_indices'];original=normalize(a[ai]);target=saved[key+'--candidate']
        u,_,vt=np.linalg.svd(original.T@(target-target.mean(0)));rotation=u@vt;translation=target.mean(0)
        assert np.max(np.abs(original@rotation+translation-target))<1e-12
        ac=a[ai].mean(0);ascale=np.sqrt(np.mean(np.sum((a[ai]-ac)**2,1)));bc=b[bi].mean(0);bscale=np.sqrt(np.mean(np.sum((b[bi]-bc)**2,1)))
        selected_a=rng.choice(len(a),8000,False);selected_b=rng.choice(len(b),8000,False)
        aa=(a[selected_a]-ac)/ascale@rotation+translation;bb=(b[selected_b]-bc)/bscale
        plots[f'candidate_{col}']=aa;plots[f'reference_{col}']=bb;plots[f'candidate_indices_{col}']=selected_a;plots[f'reference_indices_{col}']=selected_b
        info.append(dict(column=col,representation=kind,reference=name,rank=rank+1,candidate_center=ac.tolist(),candidate_scale=ascale,reference_center=bc.tolist(),reference_scale=bscale,rotation=rotation.tolist(),translation=translation.tolist(),score=row['fingerprint_distance']))
    all_xy=np.concatenate([plots[f'{kind}_{i}']@projection for i in range(3) for kind in ['candidate','reference']]);low=all_xy.min(0)-.15;high=all_xy.max(0)+.15
    fig,axes=plt.subplots(2,3,figsize=(12,8),facecolor='#071320')
    for i,item in enumerate(info):
        for j,kind in enumerate(['candidate','reference']):
            ax=axes[j,i];xy=plots[f'{kind}_{i}']@projection;ax.set_facecolor('#071320');ax.scatter(xy[:,0],xy[:,1],s=.35,c='#f6a76e' if j==0 else '#71cce5',alpha=.7,rasterized=True);ax.set_xlim(low[0],high[0]);ax.set_ylim(low[1],high[1]);ax.set_aspect('equal');ax.set_axis_off()
            label='RAAR R1051' if j==0 else item['reference'];ax.set_title(label,color='white',fontsize=12,pad=2)
        axes[0,i].text(.5,1.16,('Shape-preserving' if item['representation']=='similarity' else 'Affine-normalized')+f" / rank {item['rank']}",transform=axes[0,i].transAxes,ha='center',color='#b4c9d7',fontsize=10)
    fig.suptitle('Nearest finite-shape comparisons among 108 reference regimes',color='white',fontsize=17,y=.97)
    fig.text(.5,.055,'8,000 actual states per panel; local alignment fitted on separate 1,500-point samples. Same camera and spatial scale.',ha='center',color='#c1d2de',fontsize=9)
    fig.text(.5,.025,'Geometric resemblance and separation do not establish dynamical equivalence or historical uniqueness.',ha='center',color='#c1d2de',fontsize=9)
    fig.subplots_adjust(left=.025,right=.985,bottom=.10,top=.86,hspace=.16,wspace=.07)
    fig.savefig(out/'raar-catalog-comparison.png',dpi=230,facecolor=fig.get_facecolor());plt.close(fig)
    np.savez_compressed(out/'plotted-clouds.npz',**plots,projection=projection,low=low,high=high)
    (out/'render.json').write_text(json.dumps(dict(columns=info,points_per_panel=8000,scope='Comparison views include a different fitted orientation per pair. The affine-normalized column changes shape through covariance whitening; it is not a native-coordinate portrait.'),indent=2)+'\n',encoding='utf-8');shutil.copy2(Path(__file__),out/Path(__file__).name)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

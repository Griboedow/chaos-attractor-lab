"""Render independently simulated clouds with equal sample counts and cameras."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from compare_attractor_shapes import align

ROOT=Path(__file__).resolve().parents[1]


def main():
    base=ROOT/'results/attractor-novelty'
    summary=json.loads((base/'comparison/summary.json').read_text())
    old=np.load(ROOT/'results/twist-geometry/long-runs.npz')
    new=np.load(ROOT/'results/gaussian-vortex/long-runs.npz')
    candidates={f'Twist-{i}':z for i,z in zip(old['ids'],old['tail'])}
    candidates.update({f'Vortex-{i}':z for i,z in zip(new['ids'][1::2],new['tail'][1::2])})
    rng=np.random.default_rng(20261112)
    fig,axes=plt.subplots(3,3,figsize=(12,12),facecolor='#080c15')
    records=[];plotted={}
    for row,name in enumerate(['Vortex-77','Twist-1709','Twist-3254']):
        ref=summary['candidates'][name]['nearest'][0]['reference']
        raw=candidates[name];source=np.load(base/'reference/trajectories'/f'{ref}.npz')['points']
        a=raw[rng.choice(len(raw),8000,replace=False)];b=source[rng.choice(len(source),8000,replace=False)]
        residual,aligned=align(a[:1500],b[:1500])
        # Recover this fitted affine display transform; its linear part is a
        # uniform scale times an orthogonal matrix, with no nonlinear warping.
        transform=np.linalg.lstsq(np.column_stack((a[:1500],np.ones(1500))),aligned,rcond=None)[0]
        a=np.column_stack((a,np.ones(8000)))@transform
        bcenter=b[:1500].mean(0);brms=np.sqrt(np.mean(np.sum((b[:1500]-bcenter)**2,axis=1)))
        b=(b-bcenter)/brms
        _,basis=np.linalg.eigh(np.cov(b.T));basis=basis[:,::-1]
        az=.35;el=.25
        camera=np.array([[np.cos(az),0,np.sin(az)],[np.sin(az)*np.sin(el),np.cos(el),-np.cos(az)*np.sin(el)]])
        aa=a@basis@camera.T;bb=b@basis@camera.T
        center=(np.maximum(aa.max(0),bb.max(0))+np.minimum(aa.min(0),bb.min(0)))/2
        aa-=center;bb-=center;limit=float(max(np.abs(aa).max(),np.abs(bb).max())*1.05)
        for col,ax in enumerate(axes[row]):
            ax.set_facecolor('#080c15');ax.set_aspect('equal');ax.set_axis_off()
            ax.set_xlim(-limit,limit);ax.set_ylim(-limit,limit)
            if col in [0,2]:ax.scatter(*aa.T,s=.22,c='#ffb16b',alpha=.62,rasterized=True)
            if col in [1,2]:ax.scatter(*bb.T,s=.22,c='#62c9e9',alpha=.62,rasterized=True)
        axes[row,0].set_title(name+' / our trajectory',color='#ffb16b',fontsize=11)
        axes[row,1].set_title(ref+' / published reference',color='#62c9e9',fontsize=11)
        axes[row,2].set_title('Overlay / identical camera and scale',color='white',fontsize=10)
        records.append({'candidate':name,'reference':ref,'points_each':8000,
                        'alignment_points_each':1500,'display_fit_symmetric_rms':residual,
                        'candidate_transform':transform.tolist(),'reference_center':bcenter.tolist(),
                        'reference_rms':float(brms),'basis':basis.tolist(),'camera':camera.tolist(),
                        'plot_center':center.tolist(),'halfwidth':limit})
        plotted[name+'--candidate']=aa;plotted[name+'--reference']=bb
    fig.suptitle(f'Closest descriptor matches among {summary["reference_count"]} simulated reference regimes',color='white',fontsize=15,y=.98)
    fig.text(.5,.025,'8,000 states per panel. Rotation, reflection, translation and uniform scale allowed.\nDifferent finite clouds do not establish historical or topological novelty.',ha='center',color='#bac5d6',fontsize=11)
    fig.tight_layout(rect=(0,.055,1,.96))
    fig.savefig(ROOT/'docs/assets/attractor-novelty-comparison.png',dpi=180)
    plt.close(fig)
    (base/'comparison/render.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(base/'comparison/plotted-clouds.npz',**plotted)


if __name__=='__main__':
    main()

"""Native-state diagnostic for the certified fold-containing part of K.

The numerical centerline is a visual guide, not a certified subset of K.
The proof locates at least one fold of K inside the plotted interval hull.
"""
import argparse
import gzip
import itertools
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from propose_raar_critical_crossing import follow
from audit_raar_coverings import bounds,sha
from render_twist import camera

INK,BLUE,GOLD='#203d4b','#167b9a','#db8a24'


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.proof/'report.json').read_text());assert report['status']=='passed'
    raw=json.loads((args.proof/'input.json').read_text());p=np.array(raw['params']);S=np.array(raw['bases'][0])*np.array(raw['radii'][0]);c=np.array(raw['centers'][0])
    with gzip.open(args.proof/'volume.jsonl.gz','rt') as stream:records=[json.loads(v) for v in stream]
    critical=[v for v in records if v['classification']=='fold_if_critical'];assert critical
    low=[min(bounds(v['image'][i])[0] for v in critical) for i in range(3)]
    high=[max(bounds(v['image'][i])[1] for v in critical) for i in range(3)]
    box=np.array([[np.nextafter(float(v),-np.inf) for v in low],[np.nextafter(float(v),np.inf) for v in high]])
    u=np.linspace(-1,1,8193);starts=c+u[:,None]*S[:,0]
    saved,dets,axes=follow(starts,p,42);states=saved[-1];d=dets[-1]
    repo=Path(__file__).resolve().parents[1];refpath=repo/'results/raar-circles/render/case-1051-animation.npz'
    with np.load(refpath) as data:reference=data['states'];assert np.array_equal(data['params'],p)
    view=next(r for r in json.loads((repo/'results/raar-circles/render/analysis.json').read_text())['runs'] if r['id']==1051)
    center=np.array(view['center']);rotation=np.array(view['basis'])@camera(view['azimuth'],view['elevation']).T
    assert np.max(abs(rotation.T@rotation-np.eye(2)))<1e-12
    xy=(reference-center)@rotation;guide=(states-center)@rotation
    vertices=np.array(list(itertools.product([0,1],repeat=3)))
    corners=np.array([[box[v[i],i] for i in range(3)] for v in vertices]);outline=(corners-center)@rotation
    edges=[(i,j) for i in range(8) for j in range(i+1,8) if np.count_nonzero(vertices[i]!=vertices[j])==1]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':'#b8cbd2','xtick.color':INK,'ytick.color':INK,'pdf.fonttype':42})
    palette=LinearSegmentedColormap.from_list('context',['#ffffff','#d5e4e8','#95b9c5','#527f91'])
    fig=plt.figure(figsize=(13.5,9.0));layout=fig.add_gridspec(2,2,left=.07,right=.97,bottom=.13,top=.81,hspace=.4,height_ratios=[1.25,1])
    a,b=fig.add_subplot(layout[0,0]),fig.add_subplot(layout[0,1]);g=fig.add_subplot(layout[1,:])
    for ax,zoom in [(a,False),(b,True)]:
        lo=(outline.min(0) if zoom else xy.min(0));hi=(outline.max(0) if zoom else xy.max(0));pad=(.45 if zoom else .05)*(hi-lo)
        lo=lo-pad;hi=hi+pad
        hist,xe,ye=np.histogram2d(xy[:,0],xy[:,1],bins=(1100,850),range=list(zip(lo,hi)))
        ax.imshow(np.log1p(hist.T),origin='lower',extent=[lo[0],hi[0],lo[1],hi[1]],interpolation='nearest',cmap=palette,norm=PowerNorm(.7),alpha=.65,rasterized=True)
        for sign,color in [(-1,BLUE),(1,GOLD)]:
            take=d*sign>0
            ax.scatter(guide[take,0],guide[take,1],s=2.5 if zoom else 1.5,color=color,linewidths=0,rasterized=True)
        for i,j in edges:ax.plot(outline[[i,j],0],outline[[i,j],1],color=INK,alpha=.6,lw=.7)
        ax.set_xlim(lo[0],hi[0]);ax.set_ylim(lo[1],hi[1]);ax.set_aspect('equal');ax.axis('off')
        ax.set_title('WHOLE NUMERICAL CLOUD' if not zoom else 'BOX CONTAINING A PROVED FOLD OF K',fontsize=11,weight='bold',pad=10)
    for sign,color,label in [(-1,BLUE,'negative determinant'),(1,GOLD,'positive determinant')]:
        take=d*sign>0;g.plot(u,np.where(take,d,np.nan),color=color,lw=1.8,label=label)
    intervals=[]
    for record in critical:
        mid=np.zeros(3);rad=np.ones(3)
        for k in range(0,len(record['path']),2):
            j=int(record['path'][k]);rad[j]*=.5;mid[j]+=rad[j]*(1 if record['path'][k+1]=='1' else -1)
        intervals.append([mid[0]-rad[0],mid[0]+rad[0]])
    ulo,uhi=np.array(intervals)[:,0].min(),np.array(intervals)[:,1].max()
    g.axvspan(ulo,uhi,color=GOLD,alpha=.12,label='possible critical preimage interval')
    g.axhline(0,color=INK,lw=.7,alpha=.5);g.set_xlim(-1.02,1.02)
    g.set_xlabel('Normalized position u on the starting centerline');g.set_ylabel(r'$\det DF(F^{42}(c+uS_1))$')
    g.spines[['top','right']].set_visible(False);g.legend(loc='upper left',frameon=False,ncol=3,fontsize=9)
    fig.text(.055,.95,'A SMOOTH FOLD BELONGS TO THE ATTRACTING SET',weight='bold',fontsize=19)
    fig.text(.055,.90,'Whole-face interval signs force a critical point in K; the full window-image audit certifies its fold type.',fontsize=10.5,color='#55727f')
    fig.text(.055,.065,'Colored points: numerical images of 8,193 centerline starts after 42 steps, a guide not certified to lie in K.',fontsize=9.5,color='#55727f')
    fig.text(.055,.035,'Wire box: hull of all critical-capable image enclosures. At least one fold of K lies inside; its exact point is not located.',fontsize=9.5,color='#55727f')
    for ext in ['png','pdf']:fig.savefig(args.out/f'raar-critical-membership.{ext}',dpi=220,facecolor='white')
    plt.close(fig)
    np.savez_compressed(args.out/'figure-data.npz',params=p,normalized_u=u,starts=starts,states42=states,determinants=d,axis_minima=axes.min(0),camera_center=center,camera_rotation=rotation,critical_box=box,critical_u_intervals=np.array(intervals))
    analysis={'source_sha256':sha(Path(__file__)),'proof_report_sha256':sha(args.proof/'report.json'),'input_sha256':sha(args.proof/'input.json'),'reference_sha256':sha(refpath),'guide_points':len(u),'steps':42,'minimum_sampled_axis':float(axes.min()),'critical_capable_cells':len(critical),'critical_box_exact':[[str(v) for v in low],[str(v) for v in high]],'critical_u_hull':[float(ulo),float(uhi)],'scope':'Diagnostic, not a new attractor portrait. Native-state rigid orthographic camera. Numerical centerline not proved in K; interval hull contains at least one fold of K by the separate membership proof.'}
    (args.out/'analysis.json').write_text(json.dumps(analysis,indent=2)+'\n')
    for name in [Path(__file__).name,'propose_raar_critical_crossing.py','search_raar_circles.py','render_twist.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(analysis),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--proof',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

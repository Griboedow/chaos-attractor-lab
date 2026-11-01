"""Show all declared initial points approaching observed regimes; no connecting lines."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    cameras=json.loads(args.cameras.read_text());bg='#071320';fg='#e8f0f5';muted='#9ab0c1'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg})
    fig,axes=plt.subplots(2,5,figsize=(18,8.4),facecolor=bg);records=[]
    displayed=[0,4,16,256,4096]
    for row,ident in enumerate([1709,3254]):
        path=args.source/f'case-{ident}.npz';d=np.load(path);meta=next(c for c in cameras if c['id']==ident)
        basis=np.array(meta['pca_basis_columns']);camera=np.array(meta['projection_after_pca']);center=np.array(meta['center'])
        xy=(d['snapshots']-center)@basis@camera.T
        regular=d['return_residual'][:,2]<1e-8
        common=float(np.max(abs(xy)))*1.035;last=xy[-1]
        latecenter=(last.max(0)+last.min(0))/2;lateradius=float(np.ptp(last,axis=0).max())*.55
        bounds=[]
        for col,t in enumerate(displayed):
            ax=axes[row,col];ax.set_facecolor(bg);ax.set_aspect('equal');ax.set_axis_off()
            index=list(d['snapshot_times']).index(t);q=xy[index]
            ax.scatter(*q[~regular].T,s=.5,color='#4fd9df',alpha=.35,edgecolors='none',rasterized=True)
            if regular.any():ax.scatter(*q[regular].T,s=.5,color='#ffb65c',alpha=.35,edgecolors='none',rasterized=True)
            c=np.zeros(2) if col<4 else latecenter;r=common if col<4 else lateradius
            ax.set_xlim(c[0]-r,c[0]+r);ax.set_ylim(c[1]-r,c[1]+r)
            assert np.all(abs(q-c)<=r*(1+1e-12))
            bounds.append(dict(iterate=t,center=c.tolist(),halfwidth=r))
            title=f'n = {t:,}'
            if col==4:title+=f'   |   zoom {common/r:.1f}x'
            ax.set_title(title,fontsize=11,pad=7,color=fg)
            if col==4 and regular.any():
                cycle=np.unique(np.round(q[regular],9),axis=0)
                ax.scatter(*cycle.T,s=25,facecolors='none',edgecolors='#ffb65c',linewidths=.8)
                ax.text(.02,.02,f'3-cycle: {regular.sum():,} starts',transform=ax.transAxes,fontsize=9,color='#ffb65c')
        fig.text(.035,.83-.365*row,('Folded sail' if row==0 else 'Medusa')+f' / {ident}',fontsize=15,weight='bold')
        records.append(dict(id=ident,source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            points_per_panel=len(regular),period3_endpoint_count=int(regular.sum()),bounds=bounds,
            pca_basis_columns=basis.tolist(),projection_after_pca=camera.tolist(),center=center.tolist()))
    fig.text(.035,.95,'FROM A BALL OF STARTS TO FOLDED DYNAMICS',fontsize=22,weight='bold')
    fig.text(.035,.905,'32,768 actual states per panel. Fixed orthographic camera; common scale in the first four columns of each row.',fontsize=12,color=muted)
    fig.text(.035,.075,'Orange: final numerical three-cycle returns. Cyan: other outcomes at 4,096 iterations; this color does not certify chaos.',fontsize=11,color=muted)
    fig.text(.035,.035,'Every declared start is shown. Final panels are explicitly zoomed; no smoothing, interpolation, or connecting trajectories.',fontsize=11,color=muted)
    fig.subplots_adjust(left=.025,right=.985,top=.79,bottom=.12,wspace=.08,hspace=.38)
    fig.savefig(out/'twist-approach.png',dpi=150,facecolor=bg);plt.close(fig)
    (out/'render.json').write_text(json.dumps(dict(sources=records,classification='Each initial point is colored by its own final lag3 return rule, including earlier snapshots.',
        clipped_spatial_points=0,projection='Same rigid PCA and orthographic camera as original stills; declared camera translation and late zoom only.'),indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--cameras',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

"""Compare two motor regimes with fixed physical contrasts and a common camera."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from render_motor_four import H


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    records=[];plots=[];allbounds=[];colorbounds=[]
    for path in args.sources:
        d=np.load(path);x=d['states'];t=d['time'];p=d['parameters']
        X=np.c_[x[:,8:],-x[:,8:].sum(1)];points=X@H.T
        ids=np.flatnonzero(t>=t[-1]-1200)
        plot=points[ids];color=x[ids,:4].mean(1)
        plots.append((plot,color,p));allbounds.extend([points.min(0),points.max(0)])
        colorbounds.extend([float(color.min()),float(color.max())])
        # Compact render input contains exactly the displayed coordinates/states.
        np.savez_compressed(out/f'display-{len(records)}.npz',time=t[ids],states=x[ids],parameters=p,points=plot,color=color)
        records.append(dict(source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            source_samples=len(t),source_duration=float(t[-1]-t[0]),display_samples=len(ids),
            display_start=float(t[ids[0]]),display_end=float(t[ids[-1]]),sample_spacing=float(np.median(np.diff(t))),
            parameters=p.tolist()))
    bounds=np.array(allbounds);center=(bounds.min(0)+bounds.max(0))/2;radius=np.ptp(bounds,axis=0).max()*.55
    bg='#071320';fg='#e9f0f5';muted='#99afc0'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg})
    fig=plt.figure(figsize=(12,6.6),facecolor=bg);axes=[];artists=[]
    for i,(points,color,p) in enumerate(plots):
        ax=fig.add_axes([.015+i*.5,.12,.48,.68],projection='3d',facecolor=bg)
        seg=np.stack([points[:-1],points[1:]],axis=1)
        artist=Line3DCollection(seg,cmap='plasma',linewidth=.48,alpha=.48)
        artist.set_array((color[:-1]+color[1:])/2);artist.set_clim(min(colorbounds),max(colorbounds));ax.add_collection3d(artist)
        for setter,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):setter(c-radius,c+radius)
        ax.set_box_aspect([1,1,1],zoom=1.23);ax.set_axis_off();ax.view_init(27,-45)
        axes.append(ax);artists.append(artist)
        fig.text(.05+.5*i,.80,f"b = {p[2]:.3f}",fontsize=17,weight='bold')
    fig.text(.05,.945,'ONE MOTOR MODEL, TWO CHAOTIC GEOMETRIES',fontsize=21,weight='bold')
    fig.text(.05,.895,'Fixed displacement contrasts · common scale · eleven-dimensional dynamics',fontsize=12,color=muted)
    fig.text(.05,.075,'Final 1,200 time units of each saved history. The camera rotates around fixed data.',fontsize=11,color=muted)
    fig.text(.05,.035,'Different sampled geometries; a new attractor topology has not been established.',fontsize=10,color=muted)
    fig.savefig(out/'motor-b-regimes.png',dpi=180,facecolor=bg)
    def update(k):
        for ax in axes:ax.view_init(27,-45+360*k/72)
        return tuple(artists)
    FuncAnimation(fig,update,frames=72,interval=100,blit=False).save(out/'motor-b-regimes.gif',writer=PillowWriter(fps=10),dpi=100)
    plt.close(fig)
    (out/'render.json').write_text(json.dumps(dict(sources=records,contrast_matrix=H.tolist(),
        shared_axis_center=center.tolist(),shared_axis_radius=float(radius),bounds='Entire source histories, common bounds and equal axis units.',
        color='Mean cosine occupancy with common plasma limits.',color_limits=[min(colorbounds),max(colorbounds)],
        row_permutation='None; original labels retained in both histories.',frames=72,fps=10,
        scope='Fixed linear physical projections, no inferred 3D embedding, topology or historical novelty.'),indent=2)+'\n')
    for name in ['render_motor_b_regimes.py','render_motor_four.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,nargs=2,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

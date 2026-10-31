"""Rotate a declared final segment of the numerically checked motor ribbon."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    data=np.load(args.source);x=data['states'];t=data['time'];p=data['parameters']
    rotation=np.array([[1,-1,0],[1,1,-2],[1,1,1]],dtype=float)/np.sqrt([[2],[6],[3]])
    all_points=x[:,:3]@rotation.T
    ids=np.flatnonzero(t>=t[-1]-600)[::2];points=all_points[ids]
    lines=np.stack([points[:-1],points[1:]],axis=1)
    value=(x[ids[:-1],:3].mean(axis=1)+x[ids[1:],:3].mean(axis=1))/2
    lo=all_points.min(axis=0);hi=all_points.max(axis=0);center=(hi+lo)/2;radius=(hi-lo).max()*.55
    bg='#071320';fg='#e9f0f5';muted='#99afc0'
    fig=plt.figure(figsize=(8,8),facecolor=bg);ax=fig.add_axes([0,.12,1,.77],projection='3d',facecolor=bg)
    artist=Line3DCollection(lines,cmap='plasma',linewidth=.38,alpha=.32);artist.set_array(value);ax.add_collection3d(artist)
    for setter,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):setter(c-radius,c+radius)
    ax.set_box_aspect([1,1,1],zoom=1.15);ax.set_axis_off();ax.view_init(27,-45)
    fig.text(.055,.95,'A NARROW CHAOTIC RIBBON',fontsize=19,color=fg,weight='bold')
    fig.text(.055,.91,'Published motor model · 8D numerical candidate',fontsize=12,color=muted)
    fig.text(.055,.07,'The camera rotates; the trajectory is fixed.\nOrthogonal view of c0, c1, c2 · final 600 of 12,000 tau',fontsize=10,color=muted)
    fig.text(.055,.025,'Novelty and attractor topology are not established.',fontsize=10,color=muted)
    fig.savefig(out/'poster.png',dpi=160,facecolor=bg)
    def update(frame):
        ax.view_init(27,-45+360*frame/72)
        return (artist,)
    animation=FuncAnimation(fig,update,frames=72,interval=100,blit=False)
    animation.save(out/'rotation.gif',writer=PillowWriter(fps=10),dpi=90,savefig_kwargs={'facecolor':bg})
    plt.close(fig)
    record=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        source_samples=len(x),display_samples=len(ids),display_start=float(t[ids[0]]),display_end=float(t[ids[-1]]),
        display_rule='Final600tau, every second saved sample, uniform .04tau spacing; no smoothing or nonlinear embedding.',
        parameters=p.tolist(),rotation=rotation.tolist(),color='mean(c0,c1,c2), viridis-independent plasma color map',
        camera='72frames at10fps, elevation27deg, azimuth-45 through310deg',
        limitations='Projection of an8D flow; no topological or scientific-priority proof.')
    (out/'animation.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())

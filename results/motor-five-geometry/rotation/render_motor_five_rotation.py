"""Rotate two fixed projections of the same four-dimensional displacement set."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from render_motor_five import H


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    d=np.load(args.source);x=d['states'];t=d['time'];p=d['parameters']
    assert x.shape[1]==14
    positions=np.c_[x[:,10:],-x[:,10:].sum(1)];contrasts=positions@H.T
    keep=t>=t[-1]-800;points=contrasts[keep];color=x[keep,:5].mean(1)
    np.savez_compressed(out/'display.npz',time=t[keep],states=x[keep],parameters=p,contrasts=points,color=color)
    bg='#071320';fg='#e8f0f5';muted='#a0b3c2'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg})
    fig=plt.figure(figsize=(12,6.6),facecolor=bg);axes=[];artists=[];bounds=[]
    for panel,ids in enumerate([[0,1,2],[0,1,3]]):
        cloud=points[:,ids];allpoints=contrasts[:,ids]
        center=(allpoints.max(0)+allpoints.min(0))/2;radius=np.ptp(allpoints,axis=0).max()*.55
        bounds.append(dict(center=center.tolist(),radius=float(radius),coordinates=ids))
        ax=fig.add_axes([.015+.5*panel,.12,.48,.68],projection='3d',facecolor=bg)
        artist=Line3DCollection(np.stack([cloud[:-1],cloud[1:]],axis=1),cmap='plasma',linewidth=.5,alpha=.55)
        artist.set_array((color[:-1]+color[1:])/2);artist.set_clim(float(color.min()),float(color.max()));ax.add_collection3d(artist)
        for setter,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):setter(c-radius,c+radius)
        ax.set_box_aspect([1,1,1],zoom=1.22);ax.set_axis_off();ax.view_init(27,-45)
        axes.append(ax);artists.append(artist)
        fig.text(.05+.5*panel,.80,'cos1 / sin1 / '+('cos2' if panel==0 else 'sin2'),fontsize=15)
    fig.text(.05,.945,'FIVE MOTOR ROWS | INTERWOVEN RIBBON',fontsize=21,weight='bold')
    fig.text(.05,.895,'One trajectory, two fixed linear projections of four displacement coordinates',fontsize=12,color=muted)
    fig.text(.05,.075,'Final 800 time units of a 12,000-unit observation. The camera rotates around fixed data.',fontsize=11,color=muted)
    fig.text(.05,.035,'Fourteen-dimensional flow; neither view establishes a three-dimensional knot or a new topology.',fontsize=10,color=muted)
    fig.savefig(out/'motor-five-ribbon.png',dpi=180,facecolor=bg)
    def update(k):
        for ax in axes:ax.view_init(27,-45+360*k/72)
        return tuple(artists)
    FuncAnimation(fig,update,frames=72,interval=100,blit=False).save(out/'motor-five-ribbon.gif',writer=PillowWriter(fps=10),dpi=100)
    plt.close(fig)
    (out/'render.json').write_text(json.dumps(dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        source_samples=len(t),display_samples=int(keep.sum()),display_start=float(t[keep][0]),display_end=float(t[-1]),
        parameters=p.tolist(),contrast_matrix=H.tolist(),bounds=bounds,
        color='Mean cosine occupancy, common plasma limits in both views.',color_limits=[float(color.min()),float(color.max())],
        frames=72,fps=10,row_permutation='None',scope='Fixed projections; no 3D embedding or topology claim.'),indent=2)+'\n')
    for name in ['render_motor_five.py','render_motor_five_rotation.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

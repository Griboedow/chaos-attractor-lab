"""Fixed linear physical views and rotating four-row displacement portrait."""
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
from render_motor_four import H


def style(ax,points,elevation,azimuth,zoom=1.35):
    lo=points.min(axis=0);hi=points.max(axis=0);center=(hi+lo)/2;radius=(hi-lo).max()*.55
    for setter,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):setter(c-radius,c+radius)
    ax.set_box_aspect([1,1,1],zoom=zoom);ax.set_axis_off();ax.view_init(elevation,azimuth)


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    data=np.load(args.source);x=data['states'];t=data['time'];p=data['parameters'];ss=data['section_states']
    assert x.shape[1]==11
    positions=np.column_stack([x[:,8:],-x[:,8:].sum(axis=1)])
    points=positions@H.T;cosines=x[:,:4]@H.T;sines=x[:,4:8]@H.T
    bg='#071320';fg='#e9f0f5';muted='#99afc0'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg,'axes.labelcolor':muted,'xtick.color':muted,'ytick.color':muted})
    fig=plt.figure(figsize=(16,12),facecolor=bg)
    for row,(cloud,title,color) in enumerate([(points,'Displacements','#a2a8ff'),(cosines,'Cosine occupancies','#4bd4cd')]):
        for col,angle in enumerate([(27,-45),(27,65),(85,-90)]):
            ax=fig.add_subplot(3,3,3*row+col+1,projection='3d',facecolor=bg)
            ax.plot(*cloud[::3].T,color=color,lw=.19,alpha=.5,rasterized=True)
            style(ax,cloud,*angle);ax.set_title(f'{title} | view {col+1}',fontsize=11)
    returns=ss[:,:4]@H.T
    for col,(xx,yy,xlabel,ylabel) in enumerate([(returns[:,0],returns[:,1],'cosine contrast 1','cosine contrast 2'),
            (ss[:-1,0],ss[1:,0],'c0 at return n','c0 at return n+1'),
            (np.arange(len(data['section_times'])-1),np.diff(data['section_times']),'return index','return interval / tau')]):
        ax=fig.add_subplot(3,3,7+col,facecolor=bg);ax.scatter(xx,yy,s=2,alpha=.6,color='#f4b879')
        ax.set_xlabel(xlabel);ax.set_ylabel(ylabel);ax.spines[['top','right']].set_visible(False)
        for spine in ax.spines.values():spine.set_color('#3c5367')
        ax.set_title(['Upward X0=0 section','Scalar return plot (not a proved 1D map)','Return-time variation'][col],fontsize=10)
    fig.suptitle('FOUR MOTOR GROUPS | CAMERA AND RETURN AUDIT',x=.04,ha='left',y=.98,fontsize=19)
    fig.text(.04,.94,'Full 12,000-tau observation; every third sample in the trajectory views, all section events below.\nFixed orthonormal contrasts of four physical coordinates. 11D numerical candidate; novelty remains unresolved.',fontsize=11,color=muted,va='top')
    fig.subplots_adjust(top=.88,bottom=.07,left=.06,right=.97,hspace=.24,wspace=.27)
    fig.savefig(out/'camera-audit.png',dpi=145,facecolor=bg);plt.close(fig)
    ids=np.flatnonzero(t>=t[-1]-1200)[::2];display=points[ids]
    segments=np.stack([display[:-1],display[1:]],axis=1)
    color=(x[ids[:-1],:4].mean(axis=1)+x[ids[1:],:4].mean(axis=1))/2
    fig=plt.figure(figsize=(8,8),facecolor=bg);ax=fig.add_axes([0,.11,1,.78],projection='3d',facecolor=bg)
    artist=Line3DCollection(segments,cmap='plasma',linewidth=.42,alpha=.38);artist.set_array(color);ax.add_collection3d(artist)
    style(ax,points,27,-45,zoom=1.15)
    fig.text(.055,.95,'FOUR MOTOR GROUPS',fontsize=21,color=fg,weight='bold')
    fig.text(.055,.91,'A thin chaotic braid · 11D numerical candidate',fontsize=12,color=muted)
    fig.text(.055,.07,'The camera rotates; the trajectory is fixed.\nOrthonormal displacement contrasts · final 1,200 of 12,000 tau',fontsize=10,color=muted)
    fig.text(.055,.025,'Attractor topology and historical novelty are not established.',fontsize=10,color=muted)
    fig.savefig(out/'poster.png',dpi=160,facecolor=bg)
    def update(frame):ax.view_init(27,-45+360*frame/72);return (artist,)
    animation=FuncAnimation(fig,update,frames=72,interval=100,blit=False)
    animation.save(out/'rotation.gif',writer=PillowWriter(fps=10),dpi=100,savefig_kwargs={'facecolor':bg});plt.close(fig)
    record=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        source_samples=len(x),section_samples=len(ss),parameters=p.tolist(),contrast_matrix=H.tolist(),
        physical_displacements='X3=-(X0+X1+X2); contrast_matrix multiplies the four-vector. No nonlinear shape transformation.',
        camera_audit='Full observation, stride3; fixed cameras27/-45,27/65,85/-90; allsectionevents.',
        animation=dict(display_samples=len(ids),display_start=float(t[ids[0]]),display_end=float(t[ids[-1]]),
            rule='Final1200tau, stride2; uniform.04tau spacing. Bounds from the entire12ktau trace.',
            color='Mean of four cosine occupancies; plasma map',frames=72,fps=10,elevation=27,azimuth_start=-45),
        warning='A projection of an11D flow. Neither a proved3D embedding nor a proof of topological or historical novelty.')
    (out/'render.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    for name in ['render_motor_four_ribbon.py','render_motor_four.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

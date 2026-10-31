"""Fixed-camera and return-section audit of a slender motor candidate."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    data=np.load(args.source);states=data['states'];sections=data['section_states']
    rotation=np.array([[1/np.sqrt(2),-1/np.sqrt(2),0],[1/np.sqrt(6),1/np.sqrt(6),-2/np.sqrt(6)],
                       [1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]])
    assert np.max(abs(rotation@rotation.T-np.eye(3)))<1e-14
    bg='#071320';fg='#e2edf6';muted='#8fa8bb'
    plt.rcParams.update({'text.color':fg,'axes.labelcolor':muted,'xtick.color':muted,'ytick.color':muted})
    fig=plt.figure(figsize=(16,12),facecolor=bg)
    views=[states[:,:3]@rotation.T,states[:,3:6]@rotation.T]
    for row,points in enumerate(views):
        low=points.min(axis=0);high=points.max(axis=0);mid=(low+high)/2;radius=(high-low).max()*.53
        for col,angle in enumerate([(25,-45),(25,65),(85,-90)]):
            ax=fig.add_subplot(3,3,3*row+col+1,projection='3d',facecolor=bg)
            ax.plot(*points.T,color=['#a2a8ff','#4bd4cd'][row],lw=.18,alpha=.45,rasterized=True)
            for setter,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],mid):setter(c-radius,c+radius)
            ax.set_box_aspect([1,1,1],zoom=1.55);ax.view_init(*angle);ax.set_axis_off()
            ax.set_title(f"{'Cosine' if row==0 else 'Sine'} components | view {col+1}",fontsize=11)
    ss=sections[:,3:6]@rotation.T;cs=sections[:,:3]@rotation.T
    for col,(x,y,xlabel,ylabel) in enumerate([(cs[:,0],cs[:,1],'cosine u','cosine v'),
           (ss[:,0],ss[:,1],'sine u','sine v'),(sections[:-1,0],sections[1:,0],'c0 at return n','c0 at return n+1')]):
        ax=fig.add_subplot(3,3,7+col,facecolor=bg)
        ax.scatter(x,y,s=2,alpha=.6,color='#f4b879');ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
        ax.spines[['top','right']].set_visible(False)
        for spine in ax.spines.values():spine.set_color('#3c5367')
        ax.set_title('Upward X0=0 section' if col<2 else 'Scalar return plot (not a proved 1D map)',fontsize=10)
    fig.suptitle('SLENDER MOTOR REGIME | CAMERA AND RETURN AUDIT',x=.04,ha='left',y=.98,fontsize=19)
    fig.text(.04,.94,'All 600,000 samples / 12,000 tau. Top views use only a fixed orthogonal rotation of three source coordinates.\nA numerical 8D chaos candidate; geometric and historical novelty remain unproved.',fontsize=11,color=muted,va='top')
    fig.subplots_adjust(top=.88,bottom=.07,left=.06,right=.97,hspace=.24,wspace=.27)
    fig.savefig(out/'camera-audit.png',dpi=155,facecolor=bg);plt.close(fig)
    record=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        samples=len(states),section_samples=len(sections),rotation=rotation.tolist(),
        rotation_scope='Rigid rotation in each literal three-coordinate subspace, not a nonlinear embedding.',
        cameras=[[25,-45],[25,65],[85,-90]],section='X0=0 with positive X0 derivative, no phase filtering')
    (out/'camera-audit.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())

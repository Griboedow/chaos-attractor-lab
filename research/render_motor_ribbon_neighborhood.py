"""Shared-camera comparison of every declared local independent motor check."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.source/'analysis.json').read_text())
    assert report['declared']==report['completed']==7
    rotation=np.array([[1,-1,0],[1,1,-2],[1,1,1]],float)/np.sqrt([[2],[6],[3]])
    assert np.max(abs(rotation@rotation.T-np.eye(3)))<1e-14
    records=[];curves=[];colors=[]
    for run in report['runs']:
        path=args.source/run['file'];data=np.load(path)
        sha=hashlib.sha256(path.read_bytes()).hexdigest();assert sha==run['sha256']
        mask=data['time']>=data['time'][-1]-600
        state=data['states'][mask][::2];curves.append(state[:,:3]@rotation.T);colors.append(state[:,:3].mean(axis=1))
        records.append(dict(id=run['id'],file=path.name,sha256=sha,samples_displayed=len(state),
                            first_tau=float(data['time'][mask][0]),last_tau=float(data['time'][-1])))
    allpoints=np.vstack(curves);low=allpoints.min(axis=0);high=allpoints.max(axis=0)
    center=(low+high)/2;radius=float((high-low).max()*.53)
    bg='#081421';fg='#e8eff5';muted='#93adbf';norm=plt.Normalize(min(c.min() for c in colors),max(c.max() for c in colors))
    fig=plt.figure(figsize=(17,10),facecolor=bg)
    for i,(run,points,c) in enumerate(zip(report['runs'],curves,colors)):
        ax=fig.add_subplot(2,4,i+1,projection='3d',facecolor=bg)
        segments=np.stack([points[:-1],points[1:]],axis=1)
        lines=Line3DCollection(segments,cmap='plasma',norm=norm,linewidth=.3,alpha=.55)
        lines.set_array(c[:-1]);ax.add_collection3d(lines)
        for setter,value in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):setter(value-radius,value+radius)
        ax.set_box_aspect([1,1,1],zoom=1.28);ax.view_init(27,-45);ax.set_axis_off()
        p=run['parameters'];rate=run['spectrum'][0]
        ax.set_title(f"{'Anchor' if run['id']==0 else 'Neighbor'} {run['id']}  |  leading rate {rate:.5f}\n"
                     f"rho {p[0]:.3f}   r {p[1]:.3f}   b {p[2]:.3f}",color=fg,fontsize=10,pad=1)
    ax=fig.add_subplot(2,4,8,facecolor=bg);ax.axis('off')
    ax.text(.04,.84,'SAME VIEW, REAL COORDINATES',color=fg,fontsize=12,weight='bold')
    ax.text(.04,.71,'Fixed orthogonal rotation of c0, c1, c2.\nEqual coordinate units and shared limits.\nColor: mean cosine component.\n\nFinal 600 tau of each 3,000-tau check.\nAll seven selected checks are shown,\nincluding two weaker controls.\n\nDistinct projections do not establish\ndistinct attractor topology or priority.',color=muted,fontsize=11,linespacing=1.6,va='top')
    fig.suptitle('THE MOTOR RIBBON AND ITS NEIGHBORS',x=.045,y=.98,ha='left',color=fg,fontsize=22,weight='bold')
    fig.text(.045,.925,'Source model, eight independent states. Selection uses parameter coverage and rate controls; images were not used to choose these cases.',color=muted,fontsize=11)
    fig.subplots_adjust(left=.025,right=.985,bottom=.04,top=.87,hspace=.17,wspace=.02)
    fig.savefig(out/'gallery.png',dpi=160,facecolor=bg);plt.close(fig)
    (out/'rendering.json').write_text(json.dumps(dict(source=str(args.source),runs=records,rotation=rotation.tolist(),
        shared_center=center.tolist(),shared_radius=radius,camera=[27,-45],display='Final600tau, every second sample; no nonlinear embedding or trajectory filtering.'),indent=2)+'\n',encoding='utf8')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())

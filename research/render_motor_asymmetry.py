"""Auditable literal-coordinate gallery for the asymmetry follow-up."""
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
    records=json.loads((args.source/'analysis.json').read_text())['runs']
    bg='#071320';fg='#dbe6ef';muted='#94a8ba'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg,'axes.labelcolor':muted,'xtick.color':muted,'ytick.color':muted})
    height=4.4*len(records)+1.2
    fig=plt.figure(figsize=(15,height),facecolor=bg)
    evidence=[]
    for row,record in enumerate(records):
        path=args.source/record['file'];data=np.load(path);x=data['states'];p=data['parameters']
        views=[(x[:,3:6],['s0','s1','s2']),
               (np.column_stack([x[:,6:8],x[:,:3].mean(axis=1)]),['X0','X1','mean(c)']),
               (x[:,:3],['c0','c1','c2'])]
        for col,(points,labels) in enumerate(views):
            ax=fig.add_subplot(len(records),3,3*row+col+1,projection='3d',facecolor=bg)
            ax.plot(*points.T,color=['#43d7ce','#f4b56c','#adabff'][col],lw=.22,alpha=.47,rasterized=True)
            lo=points.min(axis=0);hi=points.max(axis=0);center=(lo+hi)/2;radius=max((hi-lo).max()/2,1e-8)*1.05
            for f,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):f(c-radius,c+radius)
            ax.set_box_aspect([1,1,1]);ax.view_init(24,-54)
            ax.set_xlabel(labels[0],labelpad=-3);ax.set_ylabel(labels[1],labelpad=-3);ax.set_zlabel(labels[2],labelpad=-3)
            ax.tick_params(labelsize=7,pad=-2);ax.grid(False)
            for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis.pane.set_edgecolor('#263748');axis.line.set_color('#394d61')
            if col==0:
                leading=', '.join(f'{v:.4f}' for v in record['spectrum'][:3])
                ax.set_title(f"Case {record['id']} | b={p[2]:g}, rho={p[0]:.4g}, r={p[1]:.4g}\nLeading rates: {leading}",fontsize=10,loc='left',pad=4)
            elif row==0:ax.set_title(['','Displacements and motor occupancy','Motor cosine components'][col],fontsize=11,pad=4)
        evidence.append(dict(file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),samples=len(x),id=record['id'],parameters=p.tolist()))
    fig.suptitle('MOTOR ASYMMETRY | PHYSICAL COORDINATE AUDIT',fontsize=19,y=1-.12/height,x=.06,ha='left')
    fig.text(.06,1-.52/height,'Three fixed linear views of each complete trace. Equal coordinate units within each panel.\n8D flow; these are projections, not reconstructed 3D attractors. No novelty claim.',fontsize=11,color=muted,va='top')
    fig.subplots_adjust(left=.025,right=.97,bottom=.025,top=1-1.3/height,hspace=.12,wspace=.06)
    fig.savefig(out/'gallery.png',dpi=160,facecolor=bg);plt.close(fig)
    (out/'gallery.json').write_text(json.dumps(dict(sources=evidence,views=['s0,s1,s2','X0,X1,mean(c)','c0,c1,c2'],
        transformation='Only listed linear observables; no nonlinear embedding, wrapping, smoothing, time selection or fitted projection.',
        camera=dict(elevation=24,azimuth=-54)),indent=2)+'\n',encoding='utf8')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())

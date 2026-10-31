"""All stratified four-row traces, fixed orthonormal physical contrasts."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

H=np.array([[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]],float)/2


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    records=json.loads((args.source/'analysis.json').read_text())['runs']
    assert all(r['N']==4 for r in records)
    assert np.max(abs(H@H.T-np.eye(3)))==0
    bg='#071320';fg='#dbe6ef';muted='#94a8ba'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg,'axes.labelcolor':muted,'xtick.color':muted,'ytick.color':muted})
    evidence=[]
    for page in range((len(records)+1)//2):
        subset=records[2*page:2*page+2];fig=plt.figure(figsize=(16,10.8),facecolor=bg)
        for row,record in enumerate(subset):
            path=args.source/record['file'];data=np.load(path);x=data['states'][::2];p=data['parameters']
            position=np.column_stack([x[:,8:],-x[:,8:].sum(axis=1)])
            views=[(position@H.T,'Displacement contrasts','#adabff'),
                   (x[:,4:8]@H.T,'Sine occupancy contrasts','#f4b56c'),
                   (x[:,:4]@H.T,'Cosine occupancy contrasts','#43d7ce')]
            for col,(points,title,color) in enumerate(views):
                ax=fig.add_subplot(2,3,3*row+col+1,projection='3d',facecolor=bg)
                ax.plot(*points.T,color=color,lw=.23,alpha=.44,rasterized=True)
                lo=points.min(axis=0);hi=points.max(axis=0);center=(lo+hi)/2;radius=max((hi-lo).max()/2,1e-8)*1.05
                for f,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):f(c-radius,c+radius)
                ax.set_box_aspect([1,1,1]);ax.view_init(27,-45);ax.grid(False)
                ax.set_xlabel('contrast 1',labelpad=-2);ax.set_ylabel('contrast 2',labelpad=-2);ax.set_zlabel('contrast 3',labelpad=-2)
                ax.tick_params(labelsize=7,pad=-2)
                for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis.pane.set_edgecolor('#263748');axis.line.set_color('#394d61')
                if row==0:ax.set_title(title,fontsize=12,pad=8)
            leading=', '.join(f'{v:.5f}' for v in record['spectrum'][:3])
            fig.text(.05,.87-.43*row,f"Case {record['id']}   |   rho {p[0]:.5g}   r {p[1]:.5g}   b {p[2]:g}   |   Leading rates {leading}",fontsize=11)
            evidence.append(dict(id=record['id'],source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                plotted_samples=len(x),page=page+1,parameters=p.tolist()))
        fig.text(.05,.965,'FOUR MOTOR ROWS | GEOMETRY SCREEN',fontsize=20)
        fig.text(.05,.93,'All selected cases, three fixed linear views; equal coordinate units in each panel.\n11-dimensional flow. These projections do not establish topological novelty.',fontsize=11,color=muted)
        fig.text(.05,.025,f'Page {page+1} of {(len(records)+1)//2}   |   Entire 1200-tau observation, every second stored sample.   |   No fitted projection, smoothing or nonlinear embedding.',fontsize=10,color=muted)
        fig.subplots_adjust(left=.02,right=.98,bottom=.055,top=.84,hspace=.25,wspace=.05)
        fig.savefig(out/f'gallery-{page+1:02}.png',dpi=140,facecolor=bg);plt.close(fig)
    (out/'gallery.json').write_text(json.dumps(dict(contrast_matrix=H.tolist(),sources=evidence,
        camera=dict(elevation=27,azimuth=-45),sampling='Entire stored observation, stride 2'),indent=2)+'\n',encoding='utf8')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

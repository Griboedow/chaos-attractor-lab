"""Fixed Fourier-coordinate gallery of every stratified five-row follow-up."""
import argparse,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

theta=2*np.pi*np.arange(5)/5
H=np.sqrt(2/5)*np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    records=json.loads((args.source/'analysis.json').read_text())['runs']
    assert all(r['N']==5 for r in records)
    assert np.max(abs(H@H.T-np.eye(4)))<1e-14 and np.max(abs(H.sum(1)))<1e-14
    bg='#071320';fg='#dbe6ef';muted='#94a8ba'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg,'axes.labelcolor':muted,'xtick.color':muted,'ytick.color':muted})
    evidence=[];pages=(len(records)+1)//2
    for page in range(pages):
        subset=records[2*page:2*page+2];fig=plt.figure(figsize=(16,10.8),facecolor=bg)
        for row,record in enumerate(subset):
            path=args.source/record['file'];data=np.load(path);x=data['states'][::2];p=data['parameters']
            position=np.column_stack([x[:,10:],-x[:,10:].sum(axis=1)]);contrast=position@H.T
            views=[(contrast[:,:3],'Displacements: cos1 / sin1 / cos2',['cos1','sin1','cos2'],'#adabff'),
                   (contrast[:,[0,1,3]],'Displacements: cos1 / sin1 / sin2',['cos1','sin1','sin2'],'#f4b56c'),
                   (x[:,5:8],'Literal sine occupancy: rows0 / 1 / 2',['s0','s1','s2'],'#43d7ce')]
            for col,(points,title,labels,color) in enumerate(views):
                ax=fig.add_subplot(2,3,3*row+col+1,projection='3d',facecolor=bg)
                ax.plot(*points.T,color=color,lw=.23,alpha=.44,rasterized=True)
                lo=points.min(axis=0);hi=points.max(axis=0);center=(lo+hi)/2;radius=max((hi-lo).max()/2,1e-8)*1.05
                for fn,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):fn(c-radius,c+radius)
                ax.set_box_aspect([1,1,1]);ax.view_init(27,-45);ax.grid(False)
                ax.set_xlabel(labels[0],labelpad=-2);ax.set_ylabel(labels[1],labelpad=-2);ax.set_zlabel(labels[2],labelpad=-2)
                ax.tick_params(labelsize=7,pad=-2)
                for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis.pane.set_edgecolor('#263748');axis.line.set_color('#394d61')
                if row==0:ax.set_title(title,fontsize=11,pad=8)
            leading=', '.join(f'{v:.5f}' for v in record['spectrum'][:3])
            fig.text(.05,.87-.43*row,f"Case {record['id']}   |   rho {p[0]:.5g}   r {p[1]:.5g}   b {p[2]:g}   |   Leading rates {leading}",fontsize=11)
            evidence.append(dict(id=record['id'],source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                plotted_samples=len(x),page=page+1,parameters=p.tolist()))
        fig.text(.05,.965,'FIVE MOTOR ROWS | GEOMETRY SCREEN',fontsize=20)
        fig.text(.05,.93,'14-dimensional flow; four independent displacement contrasts, shown in two projections.\nFixed row labels and linear observables. A 3D projection cannot preserve all four displacement coordinates.',fontsize=11,color=muted)
        fig.text(.05,.025,f'Page {page+1} of {pages}   |   Entire stored observation, every second sample.   |   Equal coordinate units within each panel; no topological novelty claim.',fontsize=10,color=muted)
        fig.subplots_adjust(left=.02,right=.98,bottom=.055,top=.84,hspace=.25,wspace=.05)
        fig.savefig(out/f'gallery-{page+1:02}.png',dpi=140,facecolor=bg);plt.close(fig)
    (out/'gallery.json').write_text(json.dumps(dict(contrast_matrix=H.tolist(),sources=evidence,
        camera=dict(elevation=27,azimuth=-45),sampling='Entire stored observation,stride2',
        limitations='Fourier coordinates are a fixed orthonormal basis of the4D zero-sum position space. Each3D panel drops a different coordinate. No reconstruction or3D group-action/topology claim.'),indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

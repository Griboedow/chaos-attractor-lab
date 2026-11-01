"""Overview, equal-unit zoom, and quantiles for the saved approach experiment."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1]
    with np.load(args.paths) as f:paths=f['states'];distance=f['relative_cloud_distances']
    with np.load(args.reference) as f:z=f['states']
    meta=json.loads(args.metadata.read_text(encoding='utf-8'));C=np.array(meta['basis'])@np.array(meta['projection']).T;center=np.array(meta['center']);xy=(z-center)@C;p=(paths-center)@C
    fig,axes=plt.subplots(1,3,figsize=(15,5.4),gridspec_kw={'width_ratios':[1,1,1.2]});norm=LogNorm(1,101)
    for ax,zoom in zip(axes[:2],[False,True]):
        ax.scatter(xy[:,0],xy[:,1],s=.1,c='#164e63',alpha=.28,rasterized=True)
        for k in range(len(paths)):ax.scatter(p[k,:101,0],p[k,:101,1],c=np.arange(101)+1,cmap='plasma_r',norm=norm,s=4 if zoom else 5,alpha=.6,edgecolors='none')
        if zoom:
            half=float(abs(xy).max()*1.3);ax.set(xlim=(-half,half),ylim=(-half,half))
        ax.set_aspect('equal');ax.set(xlabel='View coordinate 1',ylabel='View coordinate 2',title='Central zoom' if zoom else '33 initial states and their approach')
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap='plasma_r'),ax=axes[:2],orientation='horizontal',fraction=.06,pad=.14,aspect=35);cb.set_ticks([1,2,6,21,101],labels=['0','1','5','20','100']);cb.set_label('Iteration (dot color)')
    ax=axes[2];q=np.quantile(distance,[.1,.5,.9],axis=0);n=np.arange(distance.shape[1]);ax.fill_between(n,np.maximum(q[0],1e-9),q[2],color='#70a9b7',alpha=.35,label='10th–90th percentile');ax.plot(n,q[1],lw=1.1,c='#9d4939',label='Median')
    ax.set(xscale='symlog',yscale='log',xlabel='Iteration',ylabel='Distance / reference-cloud bounding-box diagonal',title='Distance to 100,000 saved reference states');ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.12)
    fig.subplots_adjust(left=.055,right=.985,top=.88,bottom=.24,wspace=.33)
    fig.text(.5,.025,'Case 988. Actual discrete states; equal spatial units in each view. The zoom clips remote approach points. Finite-cloud distance does not certify attraction.',ha='center',fontsize=8.5)
    args.figure.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.figure,dpi=160,facecolor='white');plt.close(fig)
    subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',args.figure.resolve().relative_to(root).as_posix(),'--role','diagnostic'],check=True)
    (args.out/'figure.png').write_bytes(args.figure.read_bytes());(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    report=dict(status='rendered',source_files={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.paths,args.reference,args.metadata]},figure_sha256=hashlib.sha256(args.figure.read_bytes()).hexdigest(),reference_states=len(z),approach_starts=len(paths),displayed_iterations=[0,100],scope='Sampled approach, not an invariant-set or basin certificate.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['paths','reference','metadata','figure','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())

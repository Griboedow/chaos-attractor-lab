"""Exploratory portraits of actual complete-cycle cyclic-RAAR states."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from render_portraits import COLOR


def main(args):
    args.out.parent.mkdir(parents=True,exist_ok=True);files=sorted(args.source.glob('case-*.npz'));records=[]
    fig,axes=plt.subplots(len(files),3,figsize=(13,4.3*len(files)),squeeze=False,facecolor='white')
    for row,path in enumerate(files):
        with np.load(path) as f:z=f['states'];rates=f['rates']
        center=z.mean(0);_,B=np.linalg.eigh(np.cov(z.T));B=B[:,::-1]
        if np.linalg.det(B)<0:B[:,-1]*=-1
        q=(z-center)@B;half=np.linalg.norm(q,axis=1).max()*1.04
        views=[q[:,[0,1]],q[:,[0,2]],q[:,[1,2]]]
        for col,xy in enumerate(views):
            h=np.histogram2d(xy[:,0],xy[:,1],bins=700,range=[[-half,half]]*2)[0];assert int(h.sum())==len(z)
            v=np.log1p(h);sat=np.quantile(v[h>0],.9995);im=np.minimum(v/sat,1)**.72
            ax=axes[row,col];ax.imshow(im.T,origin='lower',cmap=COLOR,vmin=0,vmax=1,interpolation='nearest');ax.set_axis_off()
            ax.set_title(f'Case {path.stem.split("-")[-1]} · view {col+1}',fontsize=12,color='#253d48')
        records.append(dict(source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),samples=len(z),center=center.tolist(),basis=B.tolist(),halfwidth=float(half),spectrum=rates.mean(0).tolist()))
    fig.suptitle('Cyclic relaxed reflections · three spatial circles',fontsize=17,color='#253d48',y=.99)
    fig.text(.5,.025,'Actual three-dimensional states after complete cycles; rigid orthogonal views and equal units.\nExploratory numerical candidates. Drawing date: 2026-11-01.',ha='center',fontsize=10,color='#253d48')
    fig.subplots_adjust(top=.93,bottom=.08,wspace=.015,hspace=.16);fig.savefig(args.out,dpi=150,facecolor='white');plt.close(fig)
    args.out.with_suffix('.json').write_text(json.dumps({'figures':records,'scope':'Exploratory shape inspection, no novelty or attractor proof.'},indent=2)+'\n',encoding='utf-8')
    root=Path(__file__).resolve().parents[1]
    if args.out.resolve().is_relative_to(root):
        subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',args.out.resolve().relative_to(root).as_posix()],check=True)
    print(args.out,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

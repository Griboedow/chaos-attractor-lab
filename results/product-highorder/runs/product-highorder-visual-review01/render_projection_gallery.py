"""Exploratory gallery with equal units and rigid cameras; every saved case shown."""
import argparse,datetime,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from render_portraits import COLOR
from render_twist import camera


def main(args):
    root=Path(__file__).resolve().parents[1];args.out.mkdir(parents=True,exist_ok=False)
    files=sorted(args.source.glob('case-*.npz'),key=lambda p:int(p.stem.split('-')[1]));records=[];projection=camera(.6,.4);date=datetime.datetime.now().astimezone().date().isoformat()
    for start in range(0,len(files),12):
        fig,axes=plt.subplots(4,3,figsize=(12,14),facecolor='white')
        for j,ax in enumerate(axes.ravel()):
            ax.set_axis_off()
            if start+j>=len(files):continue
            path=files[start+j]
            with np.load(path) as f:z=f['states'];rates=f['rates'];status=f['status'] if 'status' in f else np.array([0])
            if status[0] or not np.isfinite(z).all():ax.set_title(path.stem+' incomplete');continue
            center=z.mean(0);_,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
            if np.linalg.det(basis)<0:basis[:,-1]*=-1
            q=(z-center)@basis;half=np.linalg.norm(q,axis=1).max()*1.04;xy=q@projection.T
            h=np.histogram2d(xy[:,0],xy[:,1],bins=600,range=[[-half,half]]*2)[0];assert int(h.sum())==len(z)
            v=np.log1p(h);sat=np.quantile(v[h>0],.9995);ax.imshow(np.minimum(v/sat,1).T**.72,origin='lower',cmap=COLOR,vmin=0,vmax=1,interpolation='nearest')
            ax.set_title(f'Case {path.stem.split("-")[-1]}',fontsize=11,color='#253d48')
            records.append(dict(source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),samples=len(z),center=center.tolist(),basis=basis.tolist(),projection=projection.tolist(),halfwidth=float(half),leading_xy=float(rates[:,0].mean()),page=start//12+1))
        fig.suptitle(args.title,fontsize=16,color='#253d48',y=.98)
        fig.text(.5,.025,f'Actual R3 iterates; orthogonal camera and equal units. All saved cases, ordered by identifier.\nNumerical exploration; no novelty claim. Drawing date: {date}.',ha='center',fontsize=9,color='#253d48')
        fig.subplots_adjust(top=.94,bottom=.07,left=.025,right=.975,wspace=.02,hspace=.1);path=args.out/f'page-{start//12+1:02}.png';fig.savefig(path,dpi=130,facecolor='white');plt.close(fig)
        subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',path.resolve().relative_to(root).as_posix()],check=True)
    (args.out/'render.json').write_text(json.dumps({'records':records,'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'drawing_date':date},indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'cases':len(files),'pages':(len(files)+11)//12}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--title',required=True);main(p.parse_args())

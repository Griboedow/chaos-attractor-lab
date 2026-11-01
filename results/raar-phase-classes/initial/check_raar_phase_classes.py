"""Check whether the seven visual tips are seven separable time-residue classes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from render_twist import camera

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    with np.load(args.source) as d:z=d['states'];p=d['params']
    assert len(z)==1000000;labels=np.arange(len(z))%7;train=z[:400000];test=z[600000:];distance,nearest=cKDTree(train).query(test,workers=8);pred=labels[nearest];actual=labels[600000:];confusion=np.bincount(actual*7+pred,minlength=49).reshape(7,7);accuracy=float(np.mean(actual==pred))
    classdist=np.array([cKDTree(train[labels[:400000]==j]).query(test,workers=8)[0] for j in range(7)])
    assert np.array_equal(np.argmin(classdist,axis=0),pred)
    camrow=next(r for r in json.loads(args.camera.read_text())['runs'] if r['id']==1051);projection=np.array(camrow['basis'])@camera(camrow['azimuth'],camrow['elevation']).T;center=np.array(camrow['center']);xy=(z-center)@projection;xb=[float(xy[:,0].min()),float(xy[:,0].max())];yb=[float(xy[:,1].min()),float(xy[:,1].max())];hist=[]
    for j in range(7):
        h,_,_=np.histogram2d(xy[labels==j,0],xy[labels==j,1],bins=(700,500),range=[xb,yb]);hist.append(h.T[::-1].astype(np.uint32))
    hist=np.array(hist);assert hist.sum()==len(z);vmax=float(np.log1p(np.max(hist)));fig,axes=plt.subplots(2,4,figsize=(12,5.7));fig.suptitle('Does time modulo seven separate the arch into seven pieces?',fontsize=15)
    for j,ax in enumerate(axes.flat[:7]):ax.imshow(np.log1p(hist[j]),cmap='magma',vmin=0,vmax=vmax,extent=[*xb,*yb],aspect='equal');ax.set_title(f'Iteration index ≡ {j} (mod 7)',fontsize=10);ax.set_axis_off()
    ax=axes.flat[-1];im=ax.imshow(confusion/confusion.sum(1)[:,None],vmin=0,vmax=1,cmap='Blues');ax.set(xlabel='Nearest training phase',ylabel='Held-out phase',xticks=range(7),yticks=range(7));ax.set_title(f'3D nearest-neighbor accuracy: {accuracy:.2%}',fontsize=9);fig.colorbar(im,ax=ax,fraction=.046,pad=.04)
    fig.text(.5,.025,'Consecutive actual states; common rigid camera and density scale. First 400k train, last 400k test; 200k-step gap.',ha='center',fontsize=9);fig.subplots_adjust(left=.03,right=.98,bottom=.11,top=.89,hspace=.22,wspace=.18);fig.savefig(args.out/'raar-phase-classes.png',dpi=220);fig.savefig(args.out/'raar-phase-classes.pdf',dpi=300);plt.close(fig)
    np.savez_compressed(args.out/'phase-statistics.npz',histograms=hist,confusion=confusion,nearest_distances=distance,predicted=pred,actual=actual,class_distances=classdist,projection=projection,center=center,x_bounds=xb,y_bounds=yb)
    row=dict(source=str(args.source),source_sha256=sha(args.source),params=p.tolist(),train_indices=[0,400000],gap_indices=[400000,600000],test_indices=[600000,1000000],period=7,accuracy=accuracy,balanced_uniform_reference=1/7,confusion=confusion.tolist(),same_phase_distance_quantiles=np.quantile(classdist[actual,np.arange(len(test))],[.5,.9,.99,1]).tolist(),all_phase_distance_quantiles=np.quantile(distance,[.5,.9,.99,1]).tolist(),source_selection='One million consecutive unthinned states already frozen for R1051; no new parameter selection.',scope='Finite-sample residue-class diagnostic. Similar spatial clouds or chance-level classification do not prove ergodicity, topological mixing, absence of measurable components or asymptotic connectedness. The seven-cycle may organize local geometry without seven disjoint time-residue clouds.')
    (args.out/'report.json').write_text(json.dumps(row,indent=2),encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps(dict(accuracy=accuracy,confusion=confusion.tolist())),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=Path('results/raar-circles/render/case-1051-animation.npz'));p.add_argument('--camera',type=Path,default=Path('results/raar-circles/render/analysis.json'));p.add_argument('--out',type=Path,required=True);main(p.parse_args())

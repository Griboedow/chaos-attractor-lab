"""Diversity gallery, showing actual positive states and explicit log states."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import set_num_threads
from follow_adaptive_controller import tails
from render_twist import camera,density,colorize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    pool=np.flatnonzero(d['flags'])
    feature=[]
    for i in pool:
        z=d['tails'][i];ev=np.linalg.eigvalsh(np.cov(z.T))
        feature.append(np.r_[d['params'][i],d['rates'][i].mean(),d['counters'][i,2:]/500000,
                             np.log(np.maximum(np.std(z,axis=0),1e-10)),ev[:2]/max(ev[-1],1e-20)])
    feature=np.array(feature);feature=(feature-feature.mean(0))/np.maximum(feature.std(0),.01)
    selected=[int(np.argmax(d['rates'][pool].mean(1)))];distance=np.sum((feature-feature[selected[0]])**2,axis=1)
    while len(selected)<min(18,len(pool)):
        i=int(np.argmax(distance));assert i not in selected;selected.append(i)
        distance=np.minimum(distance,np.sum((feature-feature[i])**2,axis=1))
    chosen=pool[selected];ids=d['ids'][chosen]
    logs=tails(d['params'][chosen],d['policies'][chosen],d['ends'][chosen],200000)
    assert np.isfinite(logs).all()
    records=[]
    for j,ident in enumerate(ids):
        np.savez_compressed(args.out/f'case-{ident}.npz',log_states=logs[j],params=d['params'][chosen[j]],policy=d['policies'][chosen[j]],start=d['ends'][chosen[j]])
    for page in range((len(chosen)+2)//3):
        fig,axes=plt.subplots(3,2,figsize=(12,16),facecolor='#071320')
        for row in range(3):
            j=3*page+row
            for col in range(2):
                ax=axes[row,col];ax.set_axis_off()
                if j>=len(chosen):continue
                z=np.exp(logs[j]) if col==0 else logs[j]
                center=z.mean(0);ev,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
                if np.linalg.det(basis)<0:basis[:,-1]*=-1
                q=(z-center)@basis;shift=(q.min(0)+q.max(0))/2;q-=shift;center+=basis@shift
                half=float(max(np.linalg.norm(q,axis=1).max()*1.03,1e-12));projection=camera(.25,.3)
                h=density(q,projection,700,half);assert int(h.sum())==len(z)
                im,vmax=colorize(h);ax.imshow(im,origin='upper')
                i=chosen[j];a,b,s=d['params'][i]
                coord='(h, u, previous error)' if col==0 else 'LOGARITHMS of all three coordinates'
                ax.set_title(f'P{ids[j]} / {coord}\na={a:.4f}  b={b:.4f}  safety={s:.4f}',fontsize=10,color='#e8f0f5')
                ax.text(.025,.015,f'growth={d["rates"][i].mean():.3f}; rejection={d["counters"][i,2]/5000:.2f}%; limiter={d["counters"][i,3]/5000:.1f}%',transform=ax.transAxes,color='#b5c8d8',fontsize=9)
                records.append(dict(id=int(ids[j]),coordinates=coord,center=center.tolist(),basis=basis.tolist(),projection=projection.tolist(),halfwidth=half,vmax=vmax,points=len(z),clipped=0))
        fig.text(.04,.965,f'ADAPTIVE STEP CONTROLLER / GEOMETRY SCREEN / {page+1}',color='#e8f0f5',fontsize=18,weight='bold')
        fig.text(.04,.025,'200,000 consecutive states. Each panel uses a rigid camera and equal axis units.\nLeft: actual positive states. Right: explicitly nonlinear log-coordinate diagnostic.\nDeclared Heun/Euler controller experiment; no upstream implementation or novelty claim.',color='#a5bdce',fontsize=11)
        fig.subplots_adjust(top=.93,bottom=.085,left=.02,right=.98,wspace=.03,hspace=.15)
        fig.savefig(args.out/f'gallery-{page+1:02}.png',dpi=140,facecolor='#071320');plt.close(fig)
    report=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                eligible=len(pool),ids=ids.tolist(),selection='Largest leading exponent, then farthest standardized parameter/rate/guard-frequency/log-covariance features.',
                states='200000 next consecutive states from follow-up endpoints, all selected outcomes retained',render=records)
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    for name in ['gallery_adaptive_controller.py','follow_adaptive_controller.py','search_adaptive_controller.py','render_twist.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('Rendered',len(chosen),'cases',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

"""Reproducible diversity gallery of the tanh-Yogi scout, with rigid cameras."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from yogi_tanh_map import trajectory
from render_twist import camera,density,colorize


@njit(parallel=True)
def clouds(parameters,starts):
    states=np.empty((len(parameters),200000,3))
    for i in prange(len(parameters)):states[i]=trajectory(parameters[i],starts[i],20000,200000)
    return states


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(2)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    mask=(d['status']==0)&(d['windows'].min(1)>.005)&(np.ptp(d['tails'][:,:,0],axis=1)>1e-4)&(d['periods']==0)
    ids=np.flatnonzero(mask);features=[]
    for i in ids:
        z=d['tails'][i];cov=np.cov(z.T);e=np.maximum(np.linalg.eigvalsh(cov),1e-30);scale=np.sqrt(np.maximum(np.diag(cov),0));div=np.outer(scale,scale)
        corr=np.divide(cov,div,out=np.zeros_like(cov),where=div>1e-25)
        features.append(np.r_[d['windows'][i].mean(),np.log10(scale+1e-30),e[:2]/e[-1],corr[np.triu_indices(3,1)],np.log10(d['parameters'][i,3])])
    features=np.array(features);normed=(features-features.mean(0))/np.maximum(features.std(0),.01)
    selected=[int(np.argmax(features[:,0]))];dist=np.sum((normed-normed[selected[0]])**2,axis=1)
    while len(selected)<min(36,len(ids)):
        k=int(np.argmax(dist));assert k not in selected;selected.append(k);dist=np.minimum(dist,np.sum((normed-normed[k])**2,axis=1))
    chosen=ids[selected];p=d['parameters'][chosen];starts=d['ends'][chosen]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),selected_ids=chosen.tolist(),eligible_count=len(ids),
        selection='Largest scout leading mean first, then greedy farthest standardized leading growth, coordinate standard deviations, covariance ratios, correlations and log10(k). Both original start types remain eligible.',
        camera='Proper rigid PCA rotation, then orthographic view, equal spatial units in normalized(x,m,v); each panel fits its full extent. No nonlinear transformation.',
        warmup=20000,observed=200000,initialization='Continue from scout endpoints with autonomous map; early bias factors are numerically1 by these late times.',
        scope='Visual exploration, not fresh-start evidence. All selected clouds retained including collapse.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    z=clouds(p,starts);assert np.isfinite(z).all();np.savez_compressed(out/'histories.npz',ids=chosen,parameters=p,starts=starts,states=z)
    bg='#071320';fg='#e8f0f5';muted='#9ab0c1';records=[];pages=(len(chosen)+5)//6
    for page in range(pages):
        fig,axes=plt.subplots(2,3,figsize=(15,10.5),facecolor=bg)
        for panel,ax in enumerate(axes.ravel()):
            ax.set_axis_off();j=6*page+panel
            if j>=len(chosen):continue
            states=z[j];center=states.mean(0);ev,basis=np.linalg.eigh(np.cov(states.T));basis=basis[:,::-1]
            if np.linalg.det(basis)<0:basis[:,-1]*=-1
            q=(states-center)@basis;shift=(q.max(0)+q.min(0))/2;q-=shift;center+=basis@shift
            half=float(max(np.linalg.norm(q,axis=1).max()*1.04,1e-12));view=camera(.4)
            h=density(q,view,600,half);assert h.sum()==len(states);picture,vmax=colorize(h)
            ax.imshow(picture,origin='lower');b,c,e,k=p[j];ident=int(chosen[j])
            ax.set_title(f'S{ident} | b={b:.4f}, c={c:.4f}, e={e:.2g}, k={k:.2g}',color=fg,fontsize=10)
            ax.text(.02,.035,f'scout growth={d["windows"][ident].mean():.3f}; late x span={np.ptp(states[-4096:,0]):.3g}',transform=ax.transAxes,color=muted,fontsize=9)
            records.append(dict(id=ident,center=center.tolist(),pca_basis_columns=basis.tolist(),projection=view.tolist(),halfwidth=half,log_density_vmax=vmax,points=len(states),clipped=0))
        fig.text(.03,.95,f'TANH YOGI / EXISTING TENSORFLOW ADDONS VARIANT / {page+1} OF {pages}',color=fg,fontsize=17,weight='bold')
        fig.text(.03,.03,'200,000 actual states per panel. Rigid cameras and equal spatial units; each panel has its own extent and density scale.',color=muted,fontsize=10)
        fig.subplots_adjust(left=.02,right=.98,top=.88,bottom=.065,wspace=.025,hspace=.15)
        fig.savefig(out/f'gallery-{page+1:02}.png',dpi=140,facecolor=bg);plt.close(fig)
    (out/'render.json').write_text(json.dumps(records,indent=2)+'\n')
    for name in ['gallery_yogi_tanh.py','yogi_tanh_map.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('Rendered',len(chosen),'clouds',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

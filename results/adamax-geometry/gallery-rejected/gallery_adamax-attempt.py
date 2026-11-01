"""Inspect a declared diverse subset of the positive Adamax scout outcomes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adamax_map import trajectory
from render_twist import camera,density,colorize


@njit(parallel=True)
def histories(parameters,starts):
    result=np.empty((len(parameters),60000,3))
    for k in prange(len(parameters)):result[k]=trajectory(parameters[k],starts[k],10000,60000)
    return result


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    d=np.load(args.source);mask=(d['status']==0)&(d['windows'].min(1)>.005)&(d['extrema'][:,1,0]-d['extrema'][:,0,0]>1e-4)
    ids=np.flatnonzero(mask);features=[]
    for i in ids:
        tail=d['tails'][i];std=tail.std(0);cov=np.linalg.eigvalsh(np.cov(tail.T));cov=np.maximum(cov,1e-30)
        features.append(np.r_[d['windows'][i].mean(),d['resets'][i]/20000,np.log10(std+1e-30),cov[:2]/cov[-1],
            np.corrcoef(tail.T)[np.triu_indices(3,1)],np.mean(np.sum(np.diff(tail,axis=0)**2,axis=1))/(np.var(tail,axis=0).sum()+1e-30)])
    features=np.array(features);scale=np.std(features,axis=0);normalized=(features-features.mean(0))/np.maximum(scale,.01)
    selected=[int(np.flatnonzero(ids==i)[0]) for i in [4096,4097,4098,4099] if i in ids]
    distances=np.full(len(ids),np.inf)
    for i in selected:distances=np.minimum(distances,np.sum((normalized-normalized[i])**2,axis=1))
    while len(selected)<36:
        k=int(np.argmax(distances));selected.append(k);distances=np.minimum(distances,np.sum((normalized-normalized[k])**2,axis=1))
    selected_ids=ids[selected];parameters=d['parameters'][selected_ids];starts=d['ends'][selected_ids]
    protocol=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),selected_ids=selected_ids.tolist(),
        selection='Start with four positive declared controls; greedy farthest-point selection over standardized exponent, reset fraction, coordinate standard deviations, covariance ratios, correlations and successive-step energy.',
        eligible_count=len(ids),warmup=10000,observed=60000,camera='Rigid PCA followed by a fixed orthographic camera; equal spatial units. No nonlinear transform or separate axis scaling.',
        scope='Exploratory visual selection only; trajectories start at scout endpoints, not fresh basin samples.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    clouds=histories(parameters,starts);np.savez_compressed(out/'histories.npz',ids=selected_ids,parameters=parameters,starts=starts,states=clouds)
    bg='#071320';fg='#e8f0f5';muted='#9ab0c1';records=[]
    for page in range(6):
        fig,axes=plt.subplots(2,3,figsize=(15,10.5),facecolor=bg)
        for panel,ax in enumerate(axes.ravel()):
            k=page*6+panel;ident=int(selected_ids[k]);z=clouds[k];center=z.mean(0);ev,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
            if np.linalg.det(basis)<0:basis[:,-1]*=-1
            q=(z-center)@basis;shift=(q.max(0)+q.min(0))/2;q-=shift;center+=basis@shift
            halfwidth=float(np.max(np.linalg.norm(q,axis=1))*1.04);projection=camera(.4)
            h=density(q,projection,600,halfwidth);assert h.sum()==len(z);picture,vmax=colorize(h)
            ax.imshow(picture);ax.set_axis_off();b,c,e=parameters[k]
            ax.set_title(f'{ident} | b={b:.5f}, c={c:.5f}',color=fg,fontsize=12)
            ax.text(.02,.035,f'scout λ={d["windows"][ident].mean():.3f}, resets={d["resets"][ident]/200:.1f}%',transform=ax.transAxes,color=muted,fontsize=10)
            records.append(dict(id=ident,center=center.tolist(),pca_basis_columns=basis.tolist(),projection_after_pca=projection.tolist(),
                halfwidth=halfwidth,log_density_vmax=vmax,coordinate_extents=np.ptp(z,axis=0).tolist(),clipped_points=0))
        fig.text(.03,.95,f'ADAMAX / SCALAR QUADRATIC / GEOMETRY SCREEN {page+1}/6',color=fg,fontsize=18,weight='bold')
        fig.text(.03,.03,'60,000 actual states per panel. Rigid cameras; each panel fits its own spatial extent and log-density scale.',color=muted,fontsize=11)
        fig.subplots_adjust(left=.02,right=.98,top=.88,bottom=.065,wspace=.025,hspace=.15)
        fig.savefig(out/f'gallery-{page+1:02}.png',dpi=140,facecolor=bg);plt.close(fig)
    (out/'render.json').write_text(json.dumps(records,indent=2)+'\n')
    for name in ['gallery_adamax.py','adamax_map.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('Rendered',len(selected_ids),'candidates',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

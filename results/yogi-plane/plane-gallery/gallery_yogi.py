"""Visual follow-up of 24 diverse Yogi scout flags, including later failures."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from yogi_map import trajectory
from render_twist import camera,density,colorize


@njit(parallel=True)
def histories(parameters,starts):
    result=np.empty((len(parameters),100000,3))
    for k in prange(len(parameters)):result[k]=trajectory(parameters[k],starts[k],20000,100000)
    return result


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    archive=np.load(args.source);d={key:archive[key] for key in archive.files}
    mask=(d['status']==0)&(d['windows'].min(1)>.005)&(np.ptp(d['tails'][:,:,0],axis=1)>1e-4)&(d['periods']==0)
    ids=np.flatnonzero(mask);features=[]
    for i in ids:
        tail=d['tails'][i];std=tail.std(0);matrix=np.cov(tail.T);cov=np.maximum(np.linalg.eigvalsh(matrix),1e-30)
        denominator=np.sqrt(np.maximum(np.diag(matrix),0));denominator=np.outer(denominator,denominator)
        correlation=np.divide(matrix,denominator,out=np.zeros_like(matrix),where=denominator>1e-25)
        features.append(np.r_[d['windows'][i].mean(),d['branch_changes'][i]/20000,np.log10(std+1e-30),cov[:2]/cov[-1],
            correlation[np.triu_indices(3,1)],np.mean(np.sum(np.diff(tail,axis=0)**2,axis=1))/(np.var(tail,axis=0).sum()+1e-30)])
    features=np.array(features);assert np.isfinite(features).all();scale=np.std(features,axis=0);normalized=(features-features.mean(0))/np.maximum(scale,.01)
    selected=[int(np.argmax(features[:,0]))];distances=np.sum((normalized-normalized[selected[0]])**2,axis=1)
    while len(selected)<min(24,len(ids)):
        k=int(np.argmax(distances));assert k not in selected;selected.append(k);distances=np.minimum(distances,np.sum((normalized-normalized[k])**2,axis=1))
    selected_ids=ids[selected];parameters=d['parameters'][selected_ids];starts=d['ends'][selected_ids]
    protocol=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),selected_ids=selected_ids.tolist(),
        selection='Start with the largest scout leading mean; greedy farthest selection over standardized growth, branch-switch fraction, coordinate standard deviations, covariance ratios, correlations and successive-step energy.',
        eligible_count=len(ids),warmup=20000,observed=100000,
        camera='Proper rigid PCA rotation followed by an orthographic camera with equal spatial units in normalized (x,m,v). No nonlinear transform; v is not replaced by its square root.',
        scope='Visual selection after scout endpoints, not fresh initial-condition trials. Later collapsed or regular outcomes are retained.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    clouds=histories(parameters,starts);assert np.isfinite(clouds).all()
    np.savez_compressed(out/'histories.npz',ids=selected_ids,parameters=parameters,starts=starts,states=clouds)
    bg='#071320';fg='#e8f0f5';muted='#9ab0c1';records=[]
    for page in range((len(selected_ids)+5)//6):
        fig,axes=plt.subplots(2,3,figsize=(15,10.5),facecolor=bg)
        for panel,ax in enumerate(axes.ravel()):
            k=page*6+panel;ax.set_axis_off()
            if k>=len(selected_ids):continue
            ident=int(selected_ids[k]);z=clouds[k];center=z.mean(0);ev,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
            if np.linalg.det(basis)<0:basis[:,-1]*=-1
            q=(z-center)@basis;shift=(q.max(0)+q.min(0))/2;q-=shift;center+=basis@shift
            halfwidth=float(max(np.max(np.linalg.norm(q,axis=1))*1.04,1e-12));projection=camera(.4)
            h=density(q,projection,600,halfwidth);assert h.sum()==len(z);picture,vmax=colorize(h)
            ax.imshow(picture);b,c,e=parameters[k]
            ax.set_title(f'{ident} | b={b:.4f}, c={c:.4f}, e={e:.2g}',color=fg,fontsize=11)
            ax.text(.02,.035,f'scout growth={d["windows"][ident].mean():.3f}; tail x span={np.ptp(z[-4096:,0]):.3g}',transform=ax.transAxes,color=muted,fontsize=9)
            records.append(dict(id=ident,center=center.tolist(),pca_basis_columns=basis.tolist(),projection_after_pca=projection.tolist(),
                halfwidth=halfwidth,log_density_vmax=vmax,coordinate_extents=np.ptp(z,axis=0).tolist(),late_extents=np.ptp(z[-4096:],axis=0).tolist(),clipped_points=0))
        fig.text(.03,.95,f'YOGI / SCALAR QUADRATIC / GEOMETRY SCREEN {page+1}/4',color=fg,fontsize=18,weight='bold')
        fig.text(.03,.03,'100,000 actual states per panel. Rigid cameras; each panel fits its own extent and log-density scale.',color=muted,fontsize=11)
        fig.subplots_adjust(left=.02,right=.98,top=.88,bottom=.065,wspace=.025,hspace=.15)
        fig.savefig(out/f'gallery-{page+1:02}.png',dpi=140,facecolor=bg);plt.close(fig)
    (out/'render.json').write_text(json.dumps(records,indent=2)+'\n')
    for name in ['gallery_yogi.py','yogi_map.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('Rendered',len(selected_ids),'candidates',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

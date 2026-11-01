"""Plot every point of the declared Yogi grids and a diverse follow-up gallery."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import set_num_threads
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from gallery_yogi import histories
from render_twist import camera,density,colorize


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    with np.load(args.source/'histories.npz') as archive:d={key:archive[key] for key in archive.files}
    a=json.loads((args.source/'analysis.json').read_text());protocol=json.loads((args.source/'protocol.json').read_text());w=d['windows']
    flags=np.zeros(len(w),dtype=int);flags[np.array(a['periods'])==0]=1;flags[a['leading_ids']]=2;flags[a['two_ids']]=3
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    for s,design in enumerate(protocol['slices']):
        ids=np.flatnonzero(d['slice_ids']==s);bs=design['beta1'];cs=design['beta2'];extent=[bs[0],bs[-1],cs[0],cs[-1]]
        for j in range(3):
            ax=axes[s,j];z=(w[ids,:,j].mean(1) if j<2 else flags[ids]).reshape((32,32))
            if j==0:im=ax.imshow(z,origin='lower',extent=extent,aspect='auto',cmap='magma',vmin=0,vmax=.45,interpolation='nearest');fig.colorbar(im,ax=ax,label='Mean first branch rate',extend='both')
            elif j==1:im=ax.imshow(z,origin='lower',extent=extent,aspect='auto',cmap='RdBu_r',vmin=-.02,vmax=.02,interpolation='nearest');fig.colorbar(im,ax=ax,label='Mean second branch rate',extend='both')
            else:
                cmap=ListedColormap(['#c9cdd3','#7868ac','#dd9550','#138f91']);im=ax.imshow(z,origin='lower',extent=extent,aspect='auto',cmap=cmap,norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4),interpolation='nearest')
                cb=fig.colorbar(im,ax=ax,ticks=[0,1,2,3]);cb.ax.set_yticklabels(['Detected period','Unresolved','Leading rule','Two-direction rule'])
            ax.set(xlabel='beta1',ylabel='beta2',title=f"{design['name']} / e={design['epsilon_normalized']:.3g}")
            line=np.linspace(max(bs[0],np.sqrt(cs[0])),min(bs[-1],np.sqrt(cs[-1])),100);ax.plot(line,line*line,color='#1d5d30',ls=':',lw=1.3)
    fig.suptitle('Yogi: two cold-start parameter slices / 2,048 complete histories',fontsize=17)
    fig.supxlabel('Each point: 100,000 warmup + 200,000 observed steps. Dotted curve: beta1^2=beta2.\nRate color scales saturate at their indicated limits. Finite classification does not establish an open chaotic region.',fontsize=10)
    fig.savefig(out/'yogi-parameter-plane.png',dpi=165);plt.close(fig)
    eligible=np.array(a['leading_ids'],dtype=int);features=[]
    for i in eligible:
        z=d['tails'][i];cov=np.cov(z.T);ev=np.maximum(np.linalg.eigvalsh(cov),1e-30);sd=np.sqrt(np.maximum(np.diag(cov),0));den=np.outer(sd,sd)
        corr=np.divide(cov,den,out=np.zeros_like(cov),where=den>1e-25)
        features.append(np.r_[w[i].mean(0)[:2],np.log10(sd+1e-30),ev[:2]/ev[-1],corr[np.triu_indices(3,1)],np.mean(np.sum(np.diff(z,axis=0)**2,axis=1))/(np.var(z,axis=0).sum()+1e-30)])
    features=np.array(features);assert np.isfinite(features).all();normed=(features-features.mean(0))/np.maximum(features.std(0),.01)
    chosen=[int(np.argmax(features[:,0]))];dist=np.sum((normed-normed[chosen[0]])**2,axis=1)
    while len(chosen)<min(24,len(eligible)):
        dist[chosen]=-np.inf;k=int(np.argmax(dist));chosen.append(k);dist=np.minimum(dist,np.sum((normed-normed[k])**2,axis=1))
    ids=eligible[chosen];p=d['parameters'][ids];starts=d['ends'][ids];clouds=histories(p,starts);assert np.isfinite(clouds).all()
    np.savez_compressed(out/'histories.npz',ids=ids,parameters=p,starts=starts,states=clouds)
    records=[];bg='#071320';fg='#e8f0f5'
    for page in range((len(ids)+5)//6):
        fig,axes=plt.subplots(2,3,figsize=(15,10.5),facecolor=bg)
        for panel,ax in enumerate(axes.ravel()):
            k=page*6+panel;ax.axis('off')
            if k>=len(ids):continue
            z=clouds[k];center=z.mean(0);ev,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
            if np.linalg.det(basis)<0:basis[:,-1]*=-1
            q=(z-center)@basis;shift=(q.max(0)+q.min(0))/2;q-=shift;center+=basis@shift;radius=max(float(np.max(np.linalg.norm(q,axis=1))*1.04),1e-12)
            cam=camera(.4);h=density(q,cam,600,radius);assert h.sum()==len(z);im,vmax=colorize(h);ax.imshow(im)
            b,c,e=p[k];ax.set_title(f'Grid G{ids[k]:04d} | b={b:.4f}, c={c:.4f}, e={e:.2g}',color=fg,fontsize=11)
            ax.text(.025,.035,f'mean rates {w[ids[k],:,0].mean():+.3f}, {w[ids[k],:,1].mean():+.3f}',transform=ax.transAxes,color='#9ab0c1',fontsize=9)
            records.append(dict(grid_id=int(ids[k]),parameters=p[k].tolist(),center=center.tolist(),basis=basis.tolist(),projection=cam.tolist(),halfwidth=radius,log_density_vmax=vmax,points=len(z),clipped_points=0,late_x_span=float(np.ptp(z[-4096:,0]))))
        fig.text(.03,.95,f'YOGI / PARAMETER-GRID FOLLOW-UP / {page+1} OF 4',color=fg,fontsize=18,weight='bold')
        fig.text(.03,.03,'100,000 actual states per panel. Rigid cameras; separate extent and log-density scale per panel.',color='#9ab0c1',fontsize=11)
        fig.subplots_adjust(left=.02,right=.98,top=.88,bottom=.065,wspace=.025,hspace=.15);fig.savefig(out/f'gallery-{page+1:02d}.png',dpi=140,facecolor=bg);plt.close(fig)
    (out/'analysis.json').write_text(json.dumps(dict(source=str(args.source/'histories.npz'),source_sha256=hashlib.sha256((args.source/'histories.npz').read_bytes()).hexdigest(),
        selected_grid_ids=ids.tolist(),eligible_count=len(eligible),selection='Greedy farthest standardized feature selection, starting with largest leading mean; growth/covariance/correlation/step-energy features. No image-score or manual replacement.',
        warmup=20000,observed=100000,runs=records,scope='Exploratory visual continuation; grid identifiers are distinct from original scout case identifiers.'),indent=2)+'\n')
    for name in ['render_yogi_plane.py','gallery_yogi.py','yogi_map.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('Selected grid IDs',ids.tolist(),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

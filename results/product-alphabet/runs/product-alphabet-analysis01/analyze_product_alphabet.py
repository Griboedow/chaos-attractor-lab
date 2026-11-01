"""Time-window geometry, mean-only controls and sampled approach diagnostics."""
import argparse,hashlib,json,sys,subprocess
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from numba import njit
from search_product_alphabet import step
from analyze_soft_constellations import mean_only


@njit(cache=True)
def approaches(beta,tau,levels,order,starts,count):
    z=np.empty((len(starts),count+1,3))
    for k in range(len(starts)):
        x=starts[k].copy();z[k,0]=x
        for i in range(count):
            x,J,r,_,_,flag=step(x,beta,tau,levels,order);assert not flag;z[k,i+1]=x
    return z


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);rows=[];root=Path(__file__).resolve().parents[1]
    for name in [Path(__file__).name,'search_product_alphabet.py','analyze_soft_constellations.py','search_soft_constellations.py','search_soft_sphere.py','search_raar_constellations.py','search_spherical_constellations.py','follow_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for ident in args.ids:
        source=args.source/f'case-{ident}.npz';cohort=args.cohort/f'case-{ident}.npz';dense=args.dense/f'case-{ident}-states.npz'
        with np.load(source) as f:beta,tau=f['params'];points=f['points'];levels=f['levels'];order=int(f['order']);reference=f['states']
        with np.load(cohort) as f:starts=f['starts'];tails=f['tails'];rates=f['rates']
        with np.load(dense) as f:z=f['states']
        geometry=[]
        for k,part in enumerate(np.array_split(z,20)):
            ev=np.linalg.eigvalsh(np.cov(part.T));geometry.append(dict(window=k,mean=part.mean(0).tolist(),variance_fractions=(ev/ev.sum()).tolist(),min=part.min(0).tolist(),max=part.max(0).tolist()))
        baseline=[mean_only(x.copy(),points,tau) for x in starts]
        paths=approaches(beta,tau,levels,order,starts,1000);tree=cKDTree(reference);scale=np.linalg.norm(np.ptp(reference,axis=0));distance=tree.query(paths.reshape(-1,3),workers=1)[0].reshape(len(starts),-1)/scale
        np.savez_compressed(args.out/f'case-{ident}-approaches.npz',states=paths,starts=starts,relative_cloud_distances=distance)
        row=dict(id=ident,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cohort_sha256=hashlib.sha256(cohort.read_bytes()).hexdigest(),dense_sha256=hashlib.sha256(dense.read_bytes()).hexdigest(),dense_geometry=geometry,mean_only=[dict(initial=x.tolist(),last=y.tolist(),iterations=int(n),residual=float(r)) for x,(y,n,r) in zip(starts,baseline)],approach_distance_quantiles={str(i):np.quantile(distance[:,i],[.5,.9,1.]).tolist() for i in [0,1,2,5,10,20,50,100,200,500,1000]})
        rows.append(row)
        if ident==988:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from matplotlib.colors import LogNorm
            meta=json.loads((args.dense/f'case-{ident}.json').read_text(encoding='utf-8'));C=np.array(meta['basis'])@np.array(meta['projection']).T;center=np.array(meta['center']);xy=(reference-center)@C;pathxy=(paths-center)@C
            fig,axes=plt.subplots(1,2,figsize=(13,5.5),facecolor='white',gridspec_kw={'width_ratios':[1.1,1]})
            ax=axes[0];ax.scatter(xy[:,0],xy[:,1],s=.08,c='#164e63',alpha=.22,rasterized=True)
            for k in range(len(starts)):ax.scatter(pathxy[k,:101,0],pathxy[k,:101,1],c=np.arange(101)+1,cmap='plasma_r',norm=LogNorm(1,101),s=6,alpha=.58,edgecolors='none')
            ax.set_aspect('equal');ax.set_xlabel('Orthogonal view coordinate 1');ax.set_ylabel('Orthogonal view coordinate 2');ax.set_title('33 starts approaching the sampled cloud')
            ax=axes[1]
            for k in range(len(starts)):ax.plot(np.arange(1001),np.maximum(distance[k],1e-8),lw=.5,alpha=.25,c='#164e63')
            ax.plot(np.median(distance,axis=0),lw=1.5,c='#b65640',label='Median of 33 starts');ax.set(xscale='symlog',yscale='log',xlabel='Iteration',ylabel='Distance / reference-cloud diameter',title='Nearest sampled-state distance');ax.legend(frameon=False)
            fig.text(.5,.025,'Case 988. Dots are actual discrete states; no continuous-time interpolation. Finite-cloud distance is not a certified distance to an attractor.',ha='center',fontsize=8.5)
            fig.tight_layout(rect=[0,.05,1,.99]);args.figure.parent.mkdir(parents=True,exist_ok=True);fig.savefig(args.figure,dpi=160);plt.close(fig)
            subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',args.figure.resolve().relative_to(root).as_posix(),'--role','diagnostic'],check=True)
            row['approach_figure']=dict(path=args.figure.resolve().relative_to(root).as_posix(),sha256=hashlib.sha256(args.figure.read_bytes()).hexdigest(),samples=33*101,reference_samples=len(reference))
        print(json.dumps(dict(id=ident,mean_only_max_iterations=max(t[1] for t in baseline),first_covariance=geometry[0]['variance_fractions'],last_covariance=geometry[-1]['variance_fractions'])),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(ids=args.ids,cases=rows,scope='Finite geometry and approach diagnostics; no invariant-set identification or historical uniqueness claim.'),indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','cohort','dense','out','figure']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())

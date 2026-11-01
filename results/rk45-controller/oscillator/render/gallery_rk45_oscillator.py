"""Rigid actual-state gallery for the standard damped-oscillator RK45 map."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_rk45_oscillator import matrices,oscillator_step
from render_twist import camera,density,colorize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

@njit(cache=True)
def trajectory(params,start,warm,count):
    ms=matrices(params[0]);s=start.copy();states=np.empty((count,3))
    for t in range(warm+count):
        s,_,_,status=oscillator_step(s,params[1],*ms)
        if status:raise ValueError('Step failed')
        if t>=warm:states[t-warm]=np.array([s[0],s[1],np.exp(s[2])])
    return states

@njit(cache=True,parallel=True)
def clouds(params,starts):
    data=np.empty((len(params),200000,3))
    for i in prange(len(params)):data[i]=trajectory(params[i],starts[i],10000,200000)
    return data

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    pool=np.flatnonzero(d['flags']);features=[]
    for i in pool:
        z=d['tails'][i].copy();z[:,2]=np.exp(z[:,2]);cov=np.cov(z.T);ev=np.linalg.eigvalsh(cov);sd=np.sqrt(np.diag(cov));corr=cov/np.outer(sd,sd)
        features.append(np.r_[d['params'][i,0],np.log10(d['params'][i,1]),d['rates'][i].mean(),np.log(sd+1e-20),ev[:2]/ev[-1],corr[np.triu_indices(3,1)]])
    features=np.array(features);f=(features-features.mean(0))/np.maximum(features.std(0),.01)
    selected=[int(np.argmax(d['rates'][pool].mean(1)))];dist=np.sum((f-f[selected[0]])**2,axis=1)
    while len(selected)<min(24,len(pool)):
        i=int(np.argmax(dist));assert i not in selected;selected.append(i);dist=np.minimum(dist,np.sum((f-f[i])**2,axis=1))
    ids=np.array(args.ids,dtype=int) if args.ids else pool[selected]
    assert np.all(d['flags'][ids])
    data=clouds(d['params'][ids],d['ends'][ids]);assert np.isfinite(data).all();records=[]
    for j,i in enumerate(ids):np.savez_compressed(args.out/f'case-{i}.npz',states=data[j],params=d['params'][i],start=d['ends'][i])
    for page in range((len(ids)+5)//6):
        fig,axes=plt.subplots(2,3,figsize=(16,11),facecolor='#071320')
        for panel,ax in enumerate(axes.ravel()):
            ax.set_axis_off();j=page*6+panel
            if j>=len(ids):continue
            z=data[j];i=ids[j];center=z.mean(0);ev,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
            if np.linalg.det(basis)<0:basis[:,-1]*=-1
            q=(z-center)@basis;shift=(q.min(0)+q.max(0))/2;q-=shift;center+=basis@shift
            half=float(np.linalg.norm(q,axis=1).max()*1.04);proj=camera(.6,.4);h=density(q,proj,700,half);assert int(h.sum())==len(z)
            im,vmax=colorize(h);ax.imshow(im,origin='upper')
            ax.set_title(f'D{i} / damping={d["params"][i,0]:.6f}, rtol={d["params"][i,1]:.2g}\nScout leading rate={d["rates"][i].mean():.4f}',fontsize=10,color='#e8f0f5')
            returns=[float(np.max(np.linalg.norm(z[-2048+lag:]-z[-2048:-lag],axis=1))) for lag in range(1,513)]
            records.append(dict(id=int(i),center=center.tolist(),basis=basis.tolist(),projection=proj.tolist(),halfwidth=half,vmax=vmax,points=len(z),clipped=0,return_residuals=returns))
        fig.text(.03,.96,f'RK45 / THE DAMPED OSCILLATOR / GEOMETRY SCREEN {page+1}',fontsize=19,color='#e8f0f5',weight='bold')
        fig.text(.03,.025,'Actual (position / atol, velocity / atol, proposed step) states. 200,000 points; rigid cameras, equal spatial units.\nThe physical oscillator decays. These are numerical-controller candidates near its tolerance; no novelty claim.',fontsize=11,color='#b5c8d8')
        fig.subplots_adjust(top=.9,bottom=.08,left=.02,right=.98,wspace=.02,hspace=.16)
        fig.savefig(args.out/f'gallery-{page+1:02}.png',dpi=140,facecolor='#071320');plt.close(fig)
    report=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),warm=10000,observed=200000,eligible=len(pool),
                selection='Explicit IDs declared on the command line.' if args.ids else 'Largest leading rate followed by farthest standardized parameters, leading rate and actual-state covariance features.',render=records)
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    for n in ['gallery_rk45_oscillator.py','search_rk45_oscillator.py','render_twist.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    print('Rendered',len(ids),'cases',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+');main(p.parse_args())

"""Controlled Gaussian posterior-mean substitution in the circle reflector map.

For positive temperature this is not Euclidean projection onto a finite set.
Zero temperature is a separate, explicitly implemented hard-projection control.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_raar_constellations import step as hard_step


@njit(cache=True)
def posterior(x,points,tau):
    ds=np.sum((points-x)**2,axis=1);weights=np.exp(-(ds-ds.min())/tau);weights/=weights.sum()
    mean=np.zeros(3)
    for k in range(len(points)):mean+=weights[k]*points[k]
    covariance=np.zeros((3,3));entropy=0.
    for k in range(len(points)):
        d=points[k]-mean;covariance+=weights[k]*np.outer(d,d)
        if weights[k]>0:entropy-=weights[k]*np.log(weights[k])
    return mean,2*covariance/tau,weights.max(),entropy


@njit(cache=True)
def step(x,beta,points,tau):
    if tau==0:
        y,J,gap,r,flag=hard_step(x,beta,points)
        return y,J,r,1.,0.,flag
    mean,D,peak,entropy=posterior(x,points,tau);v=2*mean-x;r=np.hypot(v[0],v[1])
    if r<1e-12:return x,np.eye(3),r,peak,entropy,2
    u=np.array([v[0]/r,v[1]/r,0.]);t=np.array([-u[1],u[0],0.]);A=np.outer(t,t)/r
    y=beta*x+beta*u+(1-2*beta)*mean
    J=beta*np.eye(3)+beta*A@(2*D-np.eye(3))+(1-2*beta)*D
    return y,J,r,peak,entropy,0


def independent(x,beta,points,tau):
    logits=(2*points@x-np.sum(points**2,axis=1))/tau;w=np.exp(logits-logits.max());w/=w.sum();b=w@points
    rb=2*b-x;a=np.r_[rb[:2]/np.linalg.norm(rb[:2]),0.]
    return beta/2*(2*a-rb+x)+(1-beta)*b


@njit(cache=True,parallel=True)
def screen(params,points,counts,warm,observed):
    n=len(params);rates=np.full((n,10),np.nan);tails=np.full((n,2048,3),np.nan);status=np.zeros((n,2),np.int64)
    minima=np.full((n,2),np.inf);soft_stats=np.zeros((n,3));ratios=np.zeros(n)
    for k in prange(n):
        beta,tau=params[k];p=points[k,:counts[k]];x=np.array([.317,-.211,.419]);w=np.array([1.,.371,-.217]);w/=np.linalg.norm(w)
        radius=0.
        for v in p:radius=max(radius,np.linalg.norm(v))
        bound=max(np.linalg.norm(x),(beta+abs(1-2*beta)*radius)/(1-beta));sums=np.zeros(10)
        for it in range(warm+observed):
            x,J,r,peak,entropy,flag=step(x,beta,p,tau);minima[k,0]=min(minima[k,0],r);minima[k,1]=min(minima[k,1],peak)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            w=J@w;length=np.linalg.norm(w)
            if length<1e-280 or not np.isfinite(length):status[k,0]=4;break
            w/=length;ratios[k]=max(ratios[k],np.linalg.norm(x)/bound)
            if it>=warm:
                sums[(it-warm)//(observed//10)]+=np.log(length);status[k,1]+=1
                soft_stats[k,0]+=peak;soft_stats[k,1]+=entropy;soft_stats[k,2]+=peak<.99
            if it>=warm+observed-2048:tails[k,it-(warm+observed-2048)]=x
        if status[k,0]==0:rates[k]=sums/(observed//10);soft_stats[k]/=observed
    return rates,tails,status,minima,soft_stats,ratios


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    bases=[124,611,896,1504,1725,1873];relative=np.r_[0.,np.logspace(-5,0,32)];params=[];points=[];counts=[];base_ids=[];relative_tau=[];kinds=[]
    with np.load(args.source) as f:
        hard_rates=f['rates'][bases]
        for base in bases:
            p=f['points'][base,:f['counts'][base]];ds=np.sum((p[:,None,:]-p[None,:,:])**2,axis=2);np.fill_diagonal(ds,np.inf);scale=float(ds.min())
            for rel in relative:
                params.append([f['params'][base,0],rel*scale]);points.append(np.pad(p,((0,12-len(p)),(0,0))));counts.append(len(p));base_ids.append(base);relative_tau.append(rel);kinds.append(f['kinds'][base])
    params=np.array(params);points=np.array(points);counts=np.array(counts);base_ids=np.array(base_ids);relative_tau=np.array(relative_tau);kinds=np.array(kinds)
    rng=np.random.default_rng(2026110120);forms=[];diffs=[]
    for j in range(64):
        k=int(rng.integers(len(params)));p=points[k,:counts[k]];beta=params[k,0];tau=float(10**rng.uniform(-3,0));x=(p[0]+p[1])/2+rng.normal(0,.15,3)
        y,J,_,_,_,flag=step(x,beta,p,tau);assert not flag
        forms.append(float(np.max(abs(y-independent(x,beta,p,tau)))))
        h=1e-7;fd=np.column_stack([(independent(x+h*np.eye(3)[i],beta,p,tau)-independent(x-h*np.eye(3)[i],beta,p,tau))/(2*h) for i in range(3)])
        diffs.append(float(np.max(abs(J-fd))/max(1.,np.max(abs(J)))))
    assert max(forms)<1e-10 and max(diffs)<1e-5
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),bases=bases,relative_temperatures=relative.tolist(),temperature_scale='minimum squared distance between distinct codebook points',warm=5000,observed=15000,windows=10,initial=[.317,-.211,.419],map='Known RAAR algebra with finite-site projection replaced by the Gaussian posterior mean. Positive-temperature model is not standard finite-set RAAR.',soft_stats_columns=['mean maximum posterior weight','mean posterior entropy in nats','fraction with maximum posterior weight below 0.99'],checks=dict(n=64,independent_nested_formula=max(forms),relative_fd_jacobian=max(diffs)),source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Controlled regularization experiment; no new algorithm, asymptotic chaos or historical uniqueness claim.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    begin=time.time();rates,tails,status,minima,stats,ratios=screen(params,points,counts,5000,15000)
    assert np.array_equal(rates[::len(relative)],hard_rates),'Hard controls must reproduce original screen rates exactly'
    assert ratios.max()<1+1e-12;flags=(status[:,0]==0)&(rates.min(1)>.005)
    np.savez_compressed(args.out/'screen.npz',params=params,points=points,counts=counts,kinds=kinds,base_ids=base_ids,relative_tau=relative_tau,rates=rates,tails=tails,status=status,minima=minima,soft_stats=stats,boundratios=ratios,flags=flags)
    report=dict(protocol,total=len(params),complete=int((status[:,0]==0).sum()),flags=int(flags.sum()),flagged_ids=np.flatnonzero(flags).tolist(),per_base=[dict(base=base,flagged=int(flags[base_ids==base].sum())) for base in bases],hard_controls='all six spectra bitwise equal to original screen',max_boundratio=float(ratios.max()),seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

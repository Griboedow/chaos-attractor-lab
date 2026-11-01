"""Posterior-mean/sphere reflection: controlled temperature survey and follow-up.

Positive temperature is a nonprojector substitution, not finite-set RAAR.
All state coordinates are native R3; no synthetic trajectory interpolation.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_soft_constellations import posterior
from search_spherical_constellations import step as hard_step
from follow_spherical_constellations import qr_step

DEPENDENCIES=['search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py','follow_spherical_constellations.py']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n')

@njit(cache=True)
def step(x,beta,points,tau):
    if tau==0:
        y,J,gap,r,flag=hard_step(x,beta,points)
        return y,J,r,1.,0.,flag
    mean,D,peak,entropy=posterior(x,points,tau);v=2*mean-x;r=np.linalg.norm(v)
    if r<1e-12:return x,np.eye(3),r,peak,entropy,2
    u=v/r;A=(np.eye(3)-np.outer(u,u))/r
    y=beta*x+beta*u+(1-2*beta)*mean
    J=beta*np.eye(3)+beta*A@(2*D-np.eye(3))+(1-2*beta)*D
    return y,J,r,peak,entropy,0

def independent(x,beta,points,tau):
    logits=(2*points@x-np.sum(points**2,axis=1))/tau
    w=np.exp(logits-logits.max());w/=w.sum();mean=w@points
    reflected=2*mean-x;u=reflected/np.linalg.norm(reflected)
    return beta/2*(2*u-reflected+x)+(1-beta)*mean

@njit(cache=True,parallel=True)
def screen(params,points,counts,warm,observed):
    n=len(params);rates=np.full((n,10),np.nan);tails=np.full((n,2048,3),np.nan)
    status=np.zeros((n,2),np.int64);minima=np.full((n,2),np.inf);stats=np.zeros((n,3));ratios=np.zeros(n)
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
                stats[k,0]+=peak;stats[k,1]+=entropy;stats[k,2]+=peak<.99
            if it>=warm+observed-2048:tails[k,it-(warm+observed-2048)]=x
        if status[k,0]==0:rates[k]=sums/(observed//10);stats[k]/=observed
    return rates,tails,status,minima,stats,ratios

@njit(cache=True,parallel=True)
def follow(params,points,counts,warm,count):
    n=len(params);rates=np.full((n,20,3),np.nan);states=np.full((n,100000,3),np.nan)
    status=np.zeros((n,2),np.int64);minima=np.full((n,2),np.inf);stats=np.zeros((n,3))
    for k in prange(n):
        beta,tau=params[k];p=points[k,:counts[k]];x=np.array([.317,-.211,.419]);Q=np.eye(3);sums=np.zeros((20,3))
        for it in range(warm+count):
            x,J,r,peak,entropy,flag=step(x,beta,p,tau);minima[k,0]=min(minima[k,0],r);minima[k,1]=min(minima[k,1],peak)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            Q,logs,flag=qr_step(J,Q)
            if flag:status[k,0]=4;break
            if it>=warm:
                sums[(it-warm)//(count//20)]+=logs;status[k,1]+=1
                stats[k,0]+=peak;stats[k,1]+=entropy;stats[k,2]+=peak<.99
            if it>=warm+count-100000:states[k,it-(warm+count-100000)]=x
        if status[k,0]==0:rates[k]=sums/(count//20);stats[k]/=count
    return rates,states,status,minima,stats

def checks(params,points,counts):
    rng=np.random.default_rng(2026110124);forms=[];diffs=[];qr_errors=[]
    for j in range(96):
        k=int(rng.integers(len(params)));p=points[k,:counts[k]];beta=params[k,0];tau=10**rng.uniform(-3,.2);x=(p[0]+p[1])/2+rng.normal(0,.25,3)
        y,J,r,peak,entropy,flag=step(x,beta,p,tau);assert not flag
        forms.append(float(np.max(abs(y-independent(x,beta,p,tau)))))
        h=1e-7;fd=np.column_stack([(independent(x+h*np.eye(3)[i],beta,p,tau)-independent(x-h*np.eye(3)[i],beta,p,tau))/(2*h) for i in range(3)])
        diffs.append(float(np.max(abs(J-fd))/max(1.,np.max(abs(J)))))
        Q=np.linalg.qr(rng.normal(size=(3,3)))[0];W,logs,flag=qr_step(J,Q);U,R=np.linalg.qr(J@Q);U*=np.sign(np.diag(R));assert not flag
        qr_errors.append(max(float(np.max(abs(W-U))),float(np.max(abs(logs-np.log(abs(np.diag(R))))))))
    assert max(forms)<1e-10 and max(diffs)<1e-5 and max(qr_errors)<1e-9
    return dict(samples=96,independent_nested_formula=max(forms),relative_fd_jacobian=max(diffs),independent_lapack_qr=max(qr_errors))

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    for name in DEPENDENCIES:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    if args.mode=='screen':
        bases=sorted(set([124,611,896,1504,1725,1873,668,1086,1119,1346,1694,1788,1790]));relative=np.r_[0.,np.logspace(-5,0,32)]
        params=[];points=[];counts=[];kinds=[];base_ids=[];relative_tau=[]
        with np.load(args.source) as f:
            hard_rates=f['rates'][bases]
            for base in bases:
                p=f['points'][base,:f['counts'][base]];ds=np.sum((p[:,None,:]-p[None,:,:])**2,axis=2);np.fill_diagonal(ds,np.inf);scale=float(ds.min())
                for rel in relative:
                    params.append([f['params'][base,0],rel*scale]);points.append(np.pad(p,((0,12-len(p)),(0,0))));counts.append(len(p));kinds.append(f['kinds'][base]);base_ids.append(base);relative_tau.append(rel)
        params=np.array(params);points=np.array(points);counts=np.array(counts);kinds=np.array(kinds);base_ids=np.array(base_ids);relative_tau=np.array(relative_tau)
        protocol=dict(mode=args.mode,source=str(args.source),source_sha256=sha(args.source),bases=bases,base_selection='All seven hard-sphere scout flags plus all six bases of the earlier circle-temperature survey.',relative_temperatures=relative.tolist(),temperature_scale='minimum squared intersite distance',warm=5000,observed=15000,initial=[.317,-.211,.419],checks=checks(params,points,counts),scope='Controlled positive-temperature substitution in a known reflection formula; neither a new operator nor certified chaos or priority.',source_code_sha256=sha(Path(__file__)))
        save(args.out/'protocol.json',protocol);parts=[[] for _ in range(6)]
        for offset in range(0,len(params),33):
            s=slice(offset,offset+33);values=screen(params[s],points[s],counts[s],5000,15000)
            for a,b in zip(parts,values):a.append(b)
            np.savez_compressed(args.out/f'block-{offset:04}.npz',params=params[s],points=points[s],counts=counts[s],kinds=kinds[s],base_ids=base_ids[s],relative_tau=relative_tau[s],rates=values[0],tails=values[1],status=values[2],minima=values[3],soft_stats=values[4],boundratios=values[5])
            print(json.dumps(dict(completed=min(offset+33,len(params)),total=len(params),seconds=time.time()-begin)),flush=True)
        rates,tails,status,minima,stats,ratios=[np.concatenate(v) for v in parts];assert np.array_equal(rates[::len(relative)],hard_rates),'Hard controls differ from original sphere screen';assert ratios.max()<1+1e-12
        flags=(status[:,0]==0)&(rates.min(1)>.005)
        np.savez_compressed(args.out/'screen.npz',params=params,points=points,counts=counts,kinds=kinds,base_ids=base_ids,relative_tau=relative_tau,rates=rates,tails=tails,status=status,minima=minima,soft_stats=stats,boundratios=ratios,flags=flags)
        report=dict(protocol,total=len(params),complete=int((status[:,0]==0).sum()),flags=int(flags.sum()),flagged_ids=np.flatnonzero(flags).tolist(),hard_controls='All thirteen spectra bitwise equal to the original hard-sphere screen.',per_base=[dict(base=b,flagged=int(flags[base_ids==b].sum())) for b in bases],max_boundratio=float(ratios.max()),seconds=time.time()-begin)
    else:
        assert args.count>=100000 and args.count%20==0
        with np.load(args.source) as f:
            ids=np.flatnonzero(f['flags']) if args.ids is None else np.array(args.ids);params=f['params'][ids];points=f['points'][ids];counts=f['counts'][ids];kinds=f['kinds'][ids];base_ids=f['base_ids'][ids];relative_tau=f['relative_tau'][ids]
        protocol=dict(mode=args.mode,source=str(args.source),source_sha256=sha(args.source),ids=ids.tolist(),selection=args.selection or 'Every short-scout flag.',warm=args.warm,observed=args.count,initial=[.317,-.211,.419],windows=20,scope='Fresh-start finite-time numerical evidence, not an asymptotic or novelty certificate.',source_code_sha256=sha(Path(__file__)))
        save(args.out/'protocol.json',protocol);rows=[]
        for offset in range(0,len(ids),8):
            s=slice(offset,offset+8);rates,states,status,minima,stats=follow(params[s],points[s],counts[s],args.warm,args.count)
            for j,ident in enumerate(ids[s]):
                k=offset+j;np.savez_compressed(args.out/f'case-{ident}.npz',params=params[k],points=points[k,:counts[k]],kind=kinds[k],base_id=base_ids[k],relative_tau=relative_tau[k],states=states[j],rates=rates[j],status=status[j],minima=minima[j],soft_stats=stats[j])
                row=dict(id=int(ident),base_id=int(base_ids[k]),beta=float(params[k,0]),tau=float(params[k,1]),relative_tau=float(relative_tau[k]),status=status[j].tolist(),spectrum_qr_order=rates[j].mean(0).tolist(),minimum_leading_window=float(rates[j,:,0].min()),minimum_second_window=float(rates[j,:,1].min()),minimum_radius=float(minima[j,0]),soft_stats=stats[j].tolist())
                if status[j,0]==0:
                    z=states[j];ev=np.linalg.eigvalsh(np.cov(z.T));row['variance_fractions']=(ev/ev.sum()).tolist();scale=max(float(np.linalg.norm(np.ptp(z,axis=0))),1.)
                    residuals=[float(np.max(np.linalg.norm(z[-2048+lag:]-z[-2048:-lag],axis=1))) for lag in range(1,513)];periods=[k+1 for k,v in enumerate(residuals) if v<1e-8*scale];row['detected_period']=min(periods) if periods else None;row['return_residuals']=residuals
                rows.append(row)
            save(args.out/'partial.json',dict(protocol,cases=rows,seconds=time.time()-begin));print(json.dumps(dict(completed=len(rows),total=len(ids),seconds=time.time()-begin)),flush=True)
        report=dict(protocol,cases=rows,seconds=time.time()-begin)
    save(args.out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k!='cases'}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['screen','follow'],required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+');p.add_argument('--selection');p.add_argument('--count',type=int,default=500000);p.add_argument('--warm',type=int,default=30000);main(p.parse_args())

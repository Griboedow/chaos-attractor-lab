"""Declared product-alphabet posterior/reflection survey in native R3.

Equal-prior 2-, 3- and 4-level scalar alphabets factor the full 8-, 27- or
64-site Gaussian posterior. The sphere supplies the same known reflection
feedback as the canonical campaign. This is a parameter/model survey, not
a new-operator or historical-priority claim.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.stats import qmc
from follow_spherical_constellations import qr_step


@njit(cache=True)
def step(x,beta,tau,levels,count):
    mean=np.zeros(3);d=np.zeros(3);peak=1.;entropy=0.
    for j in range(3):
        logits=np.empty(count)
        for k in range(count):
            b=levels[j,k];logits[k]=(2*x[j]*b-b*b)/tau
        top=np.max(logits);weights=np.exp(logits-top);weights/=weights.sum()
        peak*=weights.max()
        for k in range(count):
            mean[j]+=weights[k]*levels[j,k]
            if weights[k]>0:entropy-=weights[k]*np.log(weights[k])
        for k in range(count):d[j]+=2/tau*weights[k]*(levels[j,k]-mean[j])**2
    v=2*mean-x;r=np.linalg.norm(v)
    if r<1e-12:return x,np.eye(3),r,peak,entropy,2
    u=v/r;y=beta*x+beta*u+(1-2*beta)*mean;J=np.empty((3,3))
    for i in range(3):
        for j in range(3):J[i,j]=beta*(i==j)+(1-2*beta)*d[i]*(i==j)+beta/r*((i==j)-u[i]*u[j])*(2*d[j]-1)
    return y,J,r,peak,entropy,0


@njit(cache=True,parallel=True)
def screen(params,levels,orders,warm,count):
    n=len(params);rates=np.full((n,10),np.nan);tails=np.full((n,2048,3),np.nan)
    status=np.zeros((n,2),np.int64);minima=np.full(n,np.inf)
    for k in prange(n):
        beta,tau=params[k];x=np.array([.317,-.211,.419]);w=np.array([1.,.371,-.217]);w/=np.linalg.norm(w);sums=np.zeros(10)
        for it in range(warm+count):
            x,J,r,_,_,flag=step(x,beta,tau,levels[k],orders[k]);minima[k]=min(minima[k],r)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            w=J@w;length=np.linalg.norm(w)
            if length<1e-280 or not np.isfinite(length):status[k,0]=4;break
            w/=length
            if it>=warm:sums[(it-warm)//(count//10)]+=np.log(length);status[k,1]+=1
            if it>=warm+count-2048:tails[k,it-(warm+count-2048)]=x
        if status[k,0]==0:rates[k]=sums/(count//10)
    return rates,tails,status,minima


@njit(cache=True,parallel=True)
def follow(params,levels,orders,warm,count):
    n=len(params);rates=np.full((n,20,3),np.nan);states=np.full((n,100000,3),np.nan)
    status=np.zeros((n,2),np.int64);minima=np.full(n,np.inf)
    for k in prange(n):
        beta,tau=params[k];x=np.array([.317,-.211,.419]);Q=np.eye(3);sums=np.zeros((20,3))
        for it in range(warm+count):
            x,J,r,_,_,flag=step(x,beta,tau,levels[k],orders[k]);minima[k]=min(minima[k],r)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            Q,logs,flag=qr_step(J,Q)
            if flag:status[k,0]=4;break
            if it>=warm:sums[(it-warm)//(count//20)]+=logs;status[k,1]+=1
            if it>=warm+count-100000:states[k,it-(warm+count-100000)]=x
        if status[k,0]==0:rates[k]=sums/(count//20)
    return rates,states,status,minima


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def checks(params,levels,orders,points,counts):
    from search_soft_sphere import independent,step as general
    rng=np.random.default_rng(2026110145);forms=[];jac=[];fd=[]
    for _ in range(192):
        k=int(rng.integers(len(params)));x=rng.uniform(-1.5,1.5,3);beta,tau=params[k];sites=points[k,:counts[k]]
        y,J,_,_,_,flag=step(x,beta,tau,levels[k],orders[k]);assert not flag
        z,K,*_=general(x,beta,sites,tau)
        forms.append(float(max(abs(y-z))));jac.append(float(np.max(abs(J-K))/max(1.,np.max(abs(K)))))
        h=1e-7;D=np.column_stack([(independent(x+h*np.eye(3)[j],beta,sites,tau)-independent(x-h*np.eye(3)[j],beta,sites,tau))/(2*h) for j in range(3)])
        fd.append(float(np.max(abs(J-D))/max(1.,np.max(abs(J)))))
    assert max(forms)<1e-10 and max(jac)<1e-9 and max(fd)<1e-5
    return dict(samples=192,max_full_posterior_formula_error=max(forms),max_full_covariance_relative_error=max(jac),max_independent_nested_finite_difference_error=max(fd))


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    for name in [Path(__file__).name,'follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    if args.mode=='screen':
        params=[];levels=[];orders=[];points=[];counts=[];rows=[]
        def add(order,beta,tau,a,c,selection,actual=None):
            grid=np.array([np.asarray(c)[j]+np.asarray(a)[j]*np.linspace(-1,1,order) for j in range(3)]) if actual is None else actual
            sites=np.array(list(itertools.product(*grid)))
            i=len(params);params.append([beta,tau]);orders.append(order);levels.append(np.pad(grid,((0,0),(0,4-order))))
            counts.append(len(sites));points.append(np.pad(sites,((0,64-len(sites)),(0,0))))
            rows.append(dict(id=i,order=order,beta=float(beta),tau=float(tau),halfwidths=np.asarray(a).tolist(),center=np.asarray(c).tolist(),selection=selection))
        with np.load(args.source) as old:
            sites=old['points'][653,:old['counts'][653]];beta,tau=old['params'][653]
            add(2,beta,tau,[.35]*3,[-.5]*3,'Exact Case 653 input control.',np.array([np.unique(sites[:,j]) for j in range(3)]))
        for order,beta,a,tau in itertools.product([2,3,4],[.9,.95,.99],[.25,.4,.65,.9],[.005,.01,.025,.06,.12,.25]):
            add(order,beta,tau,[a]*3,[0.]*3,'Complete centered, isotropic simple grid.')
        sobol=qmc.Sobol(8,scramble=True,seed=2026110144).random_base2(args.power)
        for u in sobol:
            beta=.8+.199*u[0];tau=10**(-3+2.8*u[1]);a=.15+u[2:5];c=1.6*u[5:8]-.8
            for order in [2,3,4]:add(order,beta,tau,a,c,'Same declared Sobol geometry/noise draw evaluated at every alphabet order.')
        params=np.array(params);levels=np.array(levels);orders=np.array(orders);points=np.array(points);counts=np.array(counts)
        protocol=dict(mode='screen',source_sha256=sha(args.source),seed=2026110144,sobol_power=args.power,total=len(params),warm=5000,observed=15000,initial=[.317,-.211,.419],cases=rows,checks=checks(params,levels,orders,points,counts),scope='Finite-time native-state search. Product factorization is an exact posterior identity for stored Cartesian-product inputs; no new algorithm or historical uniqueness claim.')
        save(args.out/'protocol.json',protocol);parts=[[] for _ in range(4)]
        for offset in range(0,len(params),64):
            s=slice(offset,offset+64);out=screen(params[s],levels[s],orders[s],5000,15000)
            for dest,value in zip(parts,out):dest.append(value)
            np.savez_compressed(args.out/f'block-{offset:04}.npz',params=params[s],levels=levels[s],orders=orders[s],rates=out[0],tails=out[1],status=out[2],minimum_radius=out[3])
            print(json.dumps(dict(completed=min(offset+64,len(params)),total=len(params),seconds=time.time()-begin)),flush=True)
        rates,tails,status,minima=[np.concatenate(v) for v in parts]
        positive=(status[:,0]==0)&(rates.min(1)>.005);periods=np.zeros(len(params),int)
        for i in np.flatnonzero(positive):
            z=tails[i];scale=max(float(np.linalg.norm(np.ptp(z,axis=0))),1.)
            for lag in range(1,65):
                if np.max(np.linalg.norm(z[lag:]-z[:-lag],axis=1))<1e-10*scale:periods[i]=lag;break
        flags=positive&(periods==0)
        np.savez_compressed(args.out/'screen.npz',params=params,levels=levels,orders=orders,points=points,counts=counts,kinds=orders,base_ids=np.full(len(params),-1),relative_tau=np.zeros(len(params)),rates=rates,tails=tails,status=status,minimum_radius=minima,positive_flags=positive,periods=periods,flags=flags)
        report=dict(protocol,complete=int((status[:,0]==0).sum()),positive_flags=int(positive.sum()),periodic_positive_rejections=int((positive&(periods>0)).sum()),flags=int(flags.sum()),flagged_ids=np.flatnonzero(flags).tolist(),seconds=time.time()-begin)
    else:
        assert args.count>=100000 and args.count%20==0
        with np.load(args.source) as f:
            ids=np.flatnonzero(f['flags']) if args.ids is None else np.array(args.ids)
            params=f['params'][ids];levels=f['levels'][ids];orders=f['orders'][ids];points=f['points'][ids];counts=f['counts'][ids]
        protocol=dict(mode='follow',source_sha256=sha(args.source),ids=ids.tolist(),warm=args.warm,observed=args.count,selection=args.selection or 'Every nonperiodic short-scout positive flag.',initial=[.317,-.211,.419],windows=20,scope='Fresh-start finite-time evidence, not certified chaos or priority.')
        save(args.out/'protocol.json',protocol);rows=[]
        for offset in range(0,len(ids),8):
            sl=slice(offset,offset+8);rates,states,status,minima=follow(params[sl],levels[sl],orders[sl],args.warm,args.count)
            for j,ident in enumerate(ids[sl]):
                k=offset+j
                np.savez_compressed(args.out/f'case-{ident}.npz',params=params[k],levels=levels[k],order=orders[k],points=points[k,:counts[k]],kind=orders[k],base_id=-1,relative_tau=0.,states=states[j],rates=rates[j],status=status[j],minimum_radius=minima[j])
                row=dict(id=int(ident),order=int(orders[k]),beta=float(params[k,0]),tau=float(params[k,1]),status=status[j].tolist(),minimum_radius=float(minima[j]))
                if status[j,0]==0:
                    z=states[j];ev=np.linalg.eigvalsh(np.cov(z.T));scale=max(float(np.linalg.norm(np.ptp(z,axis=0))),1.)
                    residuals=[float(np.max(np.linalg.norm(z[-2048+lag:]-z[-2048:-lag],axis=1))) for lag in range(1,513)]
                    periods=[i+1 for i,v in enumerate(residuals) if v<1e-8*scale]
                    row.update(spectrum_qr_order=rates[j].mean(0).tolist(),minimum_leading_window=float(rates[j,:,0].min()),minimum_second_window=float(rates[j,:,1].min()),variance_fractions=(ev/ev.sum()).tolist() if ev.sum()>0 else [0.,0.,0.],detected_period=min(periods) if periods else None,return_residuals=residuals)
                rows.append(row)
            save(args.out/'partial.json',dict(protocol,cases=rows,seconds=time.time()-begin));print(json.dumps(dict(completed=len(rows),total=len(ids),seconds=time.time()-begin)),flush=True)
        report=dict(protocol,cases=rows,seconds=time.time()-begin)
    save(args.out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['cases','flagged_ids']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['screen','follow'],required=True)
    for name in ['source','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--power',type=int,default=9);p.add_argument('--count',type=int,default=500000);p.add_argument('--warm',type=int,default=30000)
    p.add_argument('--ids',type=int,nargs='+');p.add_argument('--selection');main(p.parse_args())

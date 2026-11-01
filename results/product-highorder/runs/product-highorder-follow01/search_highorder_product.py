"""Declared higher-order product-posterior survey; all short flags followed."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import numpy as np
from numba import set_num_threads
from scipy.stats import qmc
from search_product_alphabet import screen,follow,step


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,r):p.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def main(a):
    a.out.mkdir(parents=True,exist_ok=False);begin=time.time();set_num_threads(4)
    names=[Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']
    for name in names:(a.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    if a.mode=='screen':
        params=[];levels=[];orders=[];rows=[]
        def add(q,beta,tau,grid,label):
            i=len(params);params.append([beta,tau]);orders.append(q);levels.append(np.pad(grid,((0,0),(0,9-q))))
            rows.append(dict(id=i,order=int(q),beta=float(beta),tau=float(tau),levels=grid.tolist(),selection=label))
        for ident in [988,489]:
            with np.load(a.source/f'case-{ident}-states.npz') as f:add(int(f['order']),*f['params'],f['levels'][:,:int(f['order'])],f'Exact earlier product Case {ident} control; local identifier differs.')
        for q,beta,half,ratio in itertools.product([5,7,9],[.90,.95,.99],[.45,.75,1.05],[.03,.15,.75]):
            grid=np.tile(half*np.linspace(-1,1,q),(3,1));tau=ratio*(2*half/(q-1))**2
            add(q,beta,tau,grid,'Complete centered isotropic grid, noise relative to spacing.')
        for u in qmc.Sobol(8,scramble=True,seed=2026110152).random_base2(8):
            beta=.85+.149*u[0];ratio=10**(-2+2.6*u[1]);half=.15+1.1*u[2:5];center=1.6*u[5:8]-.8
            for q in [5,7,9]:
                grid=center[:,None]+half[:,None]*np.linspace(-1,1,q);tau=ratio*(2*np.exp(np.log(half).mean())/(q-1))**2
                add(q,beta,tau,grid,'Same Sobol geometry and relative-noise draw at orders 5, 7 and 9.')
        params=np.array(params);levels=np.array(levels);orders=np.array(orders)
        from search_soft_sphere import step as full,independent
        rng=np.random.default_rng(2026110153);forms=[];jac=[];fd=[]
        for _ in range(192):
            i=int(rng.integers(len(params)));x=rng.uniform(-1.5,1.5,3);sites=np.array(list(itertools.product(*levels[i,:,:orders[i]])))
            y,J,*_=step(x,*params[i],levels[i],orders[i]);z,K,*_=full(x,params[i,0],sites,params[i,1]);forms.append(float(abs(y-z).max()));jac.append(float(abs(J-K).max()/max(1.,abs(J).max())))
            delta=1e-7;D=np.column_stack([(independent(x+delta*np.eye(3)[j],params[i,0],sites,params[i,1])-independent(x-delta*np.eye(3)[j],params[i,0],sites,params[i,1]))/(2*delta) for j in range(3)])
            fd.append(float(abs(J-D).max()/max(1.,abs(J).max())))
        assert max(forms)<1e-9 and max(jac)<1e-8 and max(fd)<1e-4
        protocol=dict(mode='screen',total=len(params),seed=2026110152,warm=5000,observed=15000,initial=[.317,-.211,.419],cases=rows,checks=dict(samples=192,full_site_formula_error=max(forms),full_covariance_relative_error=max(jac),nested_finite_difference_relative_error=max(fd)),scope='Exploratory visual extension to 125, 343 and 729 fixed Cartesian sites. Same established posterior/reflection operator; no new-operator or historical priority assertion.')
        save(a.out/'protocol.json',protocol);parts=[[] for _ in range(4)]
        for offset in range(0,len(params),64):
            sl=slice(offset,offset+64);v=screen(params[sl],levels[sl],orders[sl],5000,15000)
            for dst,value in zip(parts,v):dst.append(value)
            np.savez_compressed(a.out/f'block-{offset:04}.npz',params=params[sl],levels=levels[sl],orders=orders[sl],rates=v[0],tails=v[1],status=v[2],minimum_radius=v[3])
            print(json.dumps(dict(completed=min(offset+64,len(params)),total=len(params),seconds=time.time()-begin)),flush=True)
        rates,tails,status,minima=[np.concatenate(v) for v in parts];positive=(status[:,0]==0)&(rates.min(1)>.005);periods=np.zeros(len(params),int)
        for i in np.flatnonzero(positive):
            z=tails[i];scale=max(np.linalg.norm(np.ptp(z,axis=0)),1.)
            for lag in range(1,65):
                if np.max(np.linalg.norm(z[lag:]-z[:-lag],axis=1))<1e-10*scale:periods[i]=lag;break
        flags=positive&(periods==0)
        np.savez_compressed(a.out/'screen.npz',params=params,levels=levels,orders=orders,rates=rates,tails=tails,status=status,minimum_radius=minima,positive_flags=positive,periods=periods,flags=flags)
        report=dict(protocol,complete=int((status[:,0]==0).sum()),positive_flags=int(positive.sum()),periodic_positive_rejections=int((positive&(periods>0)).sum()),flags=int(flags.sum()),flagged_ids=np.flatnonzero(flags).tolist(),seconds=time.time()-begin)
    else:
        with np.load(a.source) as f:
            ids=np.flatnonzero(f['flags']) if a.ids is None else np.array(a.ids);params=f['params'][ids];levels=f['levels'][ids];orders=f['orders'][ids]
        protocol=dict(mode='follow',source_sha256=sha(a.source),ids=ids.tolist(),warm=a.warm,observed=a.count,initial=[.317,-.211,.419],selection=a.selection or 'Every nonperiodic short positive flag.',scope='Finite-time follow-up; no certified attraction or uniqueness.')
        save(a.out/'protocol.json',protocol);rows=[]
        for offset in range(0,len(ids),8):
            sl=slice(offset,offset+8);rates,states,status,minima=follow(params[sl],levels[sl],orders[sl],a.warm,a.count)
            for j,ident in enumerate(ids[sl]):
                k=offset+j;sites=np.array(list(itertools.product(*levels[k,:,:orders[k]])))
                np.savez_compressed(a.out/f'case-{ident}.npz',params=params[k],levels=levels[k],order=orders[k],points=sites,states=states[j],rates=rates[j],status=status[j],minimum_radius=minima[j])
                row=dict(id=int(ident),order=int(orders[k]),status=status[j].tolist(),minimum_radius=float(minima[j]))
                if status[j,0]==0:
                    z=states[j];ev=np.linalg.eigvalsh(np.cov(z.T));scale=max(np.linalg.norm(np.ptp(z,axis=0)),1.)
                    residuals=[float(np.max(np.linalg.norm(z[-2048+lag:]-z[-2048:-lag],axis=1))) for lag in range(1,513)]
                    periods=[i+1 for i,v in enumerate(residuals) if v<1e-8*scale]
                    row.update(spectrum=rates[j].mean(0).tolist(),minimum_leading_window=float(rates[j,:,0].min()),variance_fractions=(ev/ev.sum()).tolist() if ev.sum()>0 else [0.,0.,0.],detected_period=min(periods) if periods else None,return_residuals=residuals)
                rows.append(row)
            save(a.out/'partial.json',dict(protocol,cases=rows,seconds=time.time()-begin));print(json.dumps(dict(completed=len(rows),total=len(ids),seconds=time.time()-begin)),flush=True)
        report=dict(protocol,cases=rows,seconds=time.time()-begin)
    save(a.out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['cases','flagged_ids']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['screen','follow'],required=True)
    for name in ['source','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--count',type=int,default=500000);p.add_argument('--warm',type=int,default=30000);p.add_argument('--ids',type=int,nargs='+');p.add_argument('--selection')
    main(p.parse_args())

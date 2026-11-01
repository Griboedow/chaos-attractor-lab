"""Declared nonuniform-level visual search using the unchanged product map."""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.stats import qmc
from search_product_alphabet import screen,checks


def main(a):
    a.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    names=[Path(__file__).name,'search_highorder_product.py','search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']
    for name in names:(a.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    params=[];levels=[];orders=[];rows=[];control_sources=[]
    def add(q,p,g,label):
        i=len(params);params.append(p);levels.append(np.pad(g,((0,0),(0,9-q))));orders.append(q)
        rows.append(dict(id=i,order=q,beta=float(p[0]),tau=float(p[1]),levels=g.tolist(),selection=label))
    for old in [988,489]:
        path=a.source/f'case-{old}-states.npz'
        with np.load(path) as f:
            q=int(f['order']);add(q,f['params'],f['levels'][:,:q],f'Exact earlier product Case {old}; identifiers local to irregular survey.')
        control_sources.append(dict(case=old,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    for u in qmc.Sobol(17,scramble=True,seed=2026110154).random_base2(8):
        beta=.85+.149*u[0];ratio=10**(-2+2.6*u[1]);half=.15+1.1*u[2:5];center=1.6*u[5:8]-.8
        interior=1.8*u[8:17].reshape(3,3)-.9
        for q in [3,4,5]:
            shape=np.column_stack([-np.ones(3),np.sort(interior[:,:q-2],axis=1),np.ones(3)])
            grid=center[:,None]+half[:,None]*shape
            tau=ratio*(2*np.exp(np.log(half).mean())/(q-1))**2
            add(q,np.array([beta,tau]),grid,'Same Sobol draw at orders 3, 4 and 5; sorted nonuniform interior levels, fixed endpoints.')
    params=np.array(params);levels=np.array(levels);orders=np.array(orders);counts=orders**3;points=np.zeros((len(params),int(counts.max()),3))
    for i,q in enumerate(orders):points[i,:counts[i]]=np.array(list(itertools.product(*levels[i,:,:q])))
    protocol=dict(mode='screen',total=len(params),seed=2026110154,sobol_dimensions=17,sobol_power=8,warm=5000,observed=15000,initial=[.317,-.211,.419],control_sources=control_sources,cases=rows,checks=checks(params,levels,orders,points,counts),scope='Nonuniform quantizer levels in the same established posterior/reflection model. No new operator or historical uniqueness assertion.')
    save=lambda p,r:p.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    save(a.out/'protocol.json',protocol);parts=[[] for _ in range(4)]
    for start in range(0,len(params),64):
        sl=slice(start,start+64);values=screen(params[sl],levels[sl],orders[sl],5000,15000)
        for dest,value in zip(parts,values):dest.append(value)
        np.savez_compressed(a.out/f'block-{start:04}.npz',params=params[sl],levels=levels[sl],orders=orders[sl],rates=values[0],tails=values[1],status=values[2],minimum_radius=values[3])
    rates,tails,status,minima=[np.concatenate(v) for v in parts];positive=(status[:,0]==0)&(rates.min(1)>.005);periods=np.zeros(len(params),int)
    for i in np.flatnonzero(positive):
        z=tails[i];scale=max(1.,np.linalg.norm(np.ptp(z,axis=0)))
        for lag in range(1,65):
            if np.max(np.linalg.norm(z[lag:]-z[:-lag],axis=1))<1e-10*scale:periods[i]=lag;break
    flags=positive&(periods==0)
    np.savez_compressed(a.out/'screen.npz',params=params,levels=levels,orders=orders,points=points,counts=counts,rates=rates,tails=tails,status=status,minimum_radius=minima,positive_flags=positive,periods=periods,flags=flags)
    report=dict(protocol,complete=int((status[:,0]==0).sum()),positive_flags=int(positive.sum()),periodic_positive_rejections=int((positive&(periods>0)).sum()),flags=int(flags.sum()),flagged_ids=np.flatnonzero(flags).tolist(),seconds=time.time()-begin)
    save(a.out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['cases','flagged_ids']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())

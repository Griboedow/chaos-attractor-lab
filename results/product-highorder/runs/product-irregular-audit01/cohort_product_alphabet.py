"""Fresh-start full-spectrum and geometric controls for product alphabets."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from scipy.spatial import cKDTree
from numba import njit,prange,set_num_threads
from search_product_alphabet import step
from follow_spherical_constellations import qr_step


@njit(cache=True,parallel=True)
def cohort(beta,tau,levels,order,starts,warm,count):
    n=len(starts);rates=np.full((n,20,3),np.nan);tails=np.full((n,4096,3),np.nan)
    status=np.zeros((n,2),np.int64);minima=np.full(n,np.inf)
    for k in prange(n):
        x=starts[k].copy();Q=np.eye(3);sums=np.zeros((20,3))
        for it in range(warm+count):
            x,J,r,_,_,flag=step(x,beta,tau,levels,order);minima[k]=min(minima[k],r)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            Q,logs,flag=qr_step(J,Q)
            if flag:status[k,0]=4;break
            if it>=warm:sums[(it-warm)//(count//20)]+=logs;status[k,1]+=1
            if it>=warm+count-4096:tails[k,it-(warm+count-4096)]=x
        if status[k,0]==0:rates[k]=sums/(count//20)
    return rates,tails,status,minima


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    starts=np.vstack([[.317,-.211,.419],8*qmc.Sobol(3,scramble=True,seed=2026110146).random_base2(5)-4])
    protocol=dict(ids=args.ids,starts=starts.tolist(),seed=2026110146,warm=100000,observed=500000,windows=20,scope='Finite basin sampling; no basin, asymptotic-spectrum or uniqueness theorem.')
    def save(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    save(args.out/'protocol.json',protocol)
    for name in [Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[]
    for ident in args.ids:
        source=args.source/f'case-{ident}.npz'
        with np.load(source) as f:beta,tau=f['params'];levels=f['levels'];order=int(f['order']);reference=f['states']
        rates,tails,status,minima=cohort(beta,tau,levels,order,starts,100000,500000)
        assert (status[:,0]==0).all()
        tree=cKDTree(reference);scale=np.linalg.norm(np.ptp(reference,axis=0));details=[]
        for k,z in enumerate(tails):
            ev=np.linalg.eigvalsh(np.cov(z.T));periods=[lag for lag in range(1,513) if np.max(np.linalg.norm(z[lag:]-z[:-lag],axis=1))<1e-8*max(1.,np.linalg.norm(np.ptp(z,axis=0)))]
            details.append(dict(start=k,spectrum=rates[k].mean(0).tolist(),minimum_windows=rates[k].min(0).tolist(),variance_fractions=(ev/ev.sum()).tolist() if ev.sum()>0 else [0.,0.,0.],period=min(periods) if periods else None,reference_distance_quantiles=np.quantile(tree.query(z,workers=1)[0]/scale,[.5,.99,1.]).tolist()))
        np.savez_compressed(args.out/f'case-{ident}.npz',starts=starts,rates=rates,tails=tails,status=status,minima=minima)
        row=dict(id=ident,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),status=status.tolist(),minimum_radius=float(minima.min()),starts=details)
        rows.append(row);save(args.out/f'case-{ident}.json',row);print(json.dumps(dict(id=ident,minimum_leading=float(rates[:,:,0].min()),seconds=time.time()-begin)),flush=True)
    save(args.out/'report.json',dict(protocol,cases=rows,seconds=time.time()-begin))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())

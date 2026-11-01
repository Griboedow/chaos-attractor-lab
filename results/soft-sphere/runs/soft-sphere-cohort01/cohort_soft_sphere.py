"""Fresh-start basin sampling for the selected posterior-mean sphere maps."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from scipy.spatial import cKDTree
from numba import njit,prange,set_num_threads
from search_soft_sphere import step


@njit(cache=True,parallel=True)
def cohort(beta,tau,points,starts,warm,count):
    n=len(starts);rates=np.full((n,20),np.nan);tails=np.full((n,4096,3),np.nan)
    status=np.zeros((n,2),np.int64);minima=np.full(n,np.inf)
    for k in prange(n):
        x=starts[k].copy();w=np.array([1.,.371,-.217]);w/=np.linalg.norm(w);sums=np.zeros(20)
        for it in range(warm+count):
            x,j,r,peak,entropy,flag=step(x,beta,points,tau);minima[k]=min(minima[k],r)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            w=j@w;length=np.linalg.norm(w)
            if length<1e-280 or not np.isfinite(length):status[k,0]=4;break
            w/=length
            if it>=warm:sums[(it-warm)//(count//20)]+=np.log(length);status[k,1]+=1
            if it>=warm+count-4096:tails[k,it-(warm+count-4096)]=x
        if status[k,0]==0:rates[k]=sums/(count//20)
    return rates,tails,status,minima


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    starts=np.vstack([[.317,-.211,.419],8*qmc.Sobol(3,scramble=True,seed=2026110121).random_base2(5)-4])
    protocol=dict(ids=args.ids,starts=starts.tolist(),warm=100000,observed=500000,windows=20,purpose='Finite fresh-start robustness; no certified basin or unique invariant measure.',distance='Last 4096 states to the reference 100000-state cloud, divided by its bounding-box diagonal; sampling diagnostic only.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py','follow_spherical_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[];begin=time.time()
    for ident in args.ids:
        source=args.source/f'case-{ident}.npz'
        with np.load(source) as f:beta,tau=f['params'];points=f['points'];reference=f['states']
        rates,tails,status,minima=cohort(beta,tau,points,starts,100000,500000)
        distances=np.full((len(starts),3),np.nan);tree=cKDTree(reference);scale=np.linalg.norm(np.ptp(reference,axis=0))
        for k in range(len(starts)):
            if status[k,0]==0:distances[k]=np.quantile(tree.query(tails[k],workers=1)[0]/scale,[.5,.99,1.])
        np.savez_compressed(args.out/f'case-{ident}.npz',starts=starts,rates=rates,tails=tails,status=status,minima=minima,distances=distances)
        row=dict(id=ident,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),status=status.tolist(),mean_rates=rates.mean(1).tolist(),minimum_windows=rates.min(1).tolist(),minimum_radii=minima.tolist(),distance_quantiles=distances.tolist())
        rows.append(row);(args.out/f'case-{ident}.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(json.dumps(dict(id=ident,complete=int((status[:,0]==0).sum()),minimum_window=float(np.nanmin(rates)),seconds=time.time()-begin)),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(protocol,cases=rows,seconds=time.time()-begin),indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())

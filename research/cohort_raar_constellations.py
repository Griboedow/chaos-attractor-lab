"""Declared fresh-start cohort for selected circle / finite-point RAAR maps."""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from scipy.spatial import cKDTree
from numba import njit, prange, set_num_threads
from search_raar_constellations import step


@njit(cache=True, parallel=True)
def cohort(beta, points, starts, warm, count):
    n=len(starts); rates=np.full((n,20),np.nan); tails=np.full((n,4096,3),np.nan)
    status=np.zeros((n,2),np.int64); minima=np.full((n,2),np.inf)
    visits=np.zeros((n,len(points)),np.int64)
    for k in prange(n):
        x=starts[k].copy(); w=np.array([1.,.371,-.217]); w/=np.linalg.norm(w); sums=np.zeros(20)
        for it in range(warm+count):
            branch=np.argmin(np.sum((points-x)**2,axis=1))
            x,j,gap,axis,flag=step(x,beta,points)
            minima[k,0]=min(minima[k,0],gap); minima[k,1]=min(minima[k,1],axis)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            w=j@w; length=np.linalg.norm(w)
            if length<1e-280 or not np.isfinite(length):status[k,0]=4;break
            w/=length
            if it>=warm:
                sums[(it-warm)//(count//20)]+=np.log(length);status[k,1]+=1;visits[k,branch]+=1
            if it>=warm+count-4096:tails[k,it-(warm+count-4096)]=x
        if status[k,0]==0:rates[k]=sums/(count//20)
    return rates,tails,status,minima,visits


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6)
    starts=np.vstack([[.317,-.211,.419],8*qmc.Sobol(3,scramble=True,seed=2026110117).random_base2(5)-4])
    protocol=dict(ids=args.ids,starts=starts.tolist(),warm=100000,observed=500000,windows=20,
                  purpose='Fresh-start robustness; finite samples do not prove a basin or a unique invariant measure.',
                  distance='Last 4096 states to the reference 100000-state cloud, normalized by its bounding-box diagonal; sampling diagnostic only.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    for name in [Path(__file__).name,'search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[];begin=time.time()
    for ident in args.ids:
        source=args.source/f'case-{ident}.npz'
        with np.load(source) as f:beta=float(f['params'][0]);points=f['points'];reference=f['states']
        rates,tails,status,minima,visits=cohort(beta,points,starts,100000,500000)
        distances=np.full((len(starts),3),np.nan);tree=cKDTree(reference);scale=np.linalg.norm(np.ptp(reference,axis=0))
        for k in range(len(starts)):
            if status[k,0]==0:distances[k]=np.quantile(tree.query(tails[k],workers=1)[0]/scale,[.5,.99,1.])
        np.savez_compressed(args.out/f'case-{ident}.npz',starts=starts,rates=rates,tails=tails,status=status,minima=minima,visits=visits,distances=distances)
        row=dict(id=ident,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),status=status.tolist(),mean_rates=rates.mean(1).tolist(),minimum_windows=rates.min(1).tolist(),minima=minima.tolist(),visits=visits.tolist(),distance_quantiles=distances.tolist())
        rows.append(row);(args.out/f'case-{ident}.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(dict(id=ident,complete=int((status[:,0]==0).sum()),minimum_window=float(np.nanmin(rates)),seconds=time.time()-begin)),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(protocol,cases=rows,seconds=time.time()-begin),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())

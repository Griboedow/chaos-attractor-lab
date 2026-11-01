"""Full local spectra and fresh-start clouds for every soft-scout flag."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_soft_constellations import step
from follow_spherical_constellations import qr_step


@njit(cache=True,parallel=True)
def follow(params,points,counts,warm,count):
    n=len(params);rates=np.full((n,20,3),np.nan);states=np.full((n,100000,3),np.nan)
    status=np.zeros((n,2),np.int64);minima=np.full((n,2),np.inf);soft_stats=np.zeros((n,3))
    for k in prange(n):
        beta,tau=params[k];p=points[k,:counts[k]];x=np.array([.317,-.211,.419]);Q=np.eye(3);sums=np.zeros((20,3))
        for it in range(warm+count):
            x,J,r,peak,entropy,flag=step(x,beta,p,tau);minima[k,0]=min(minima[k,0],r);minima[k,1]=min(minima[k,1],peak)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            Q,logs,flag=qr_step(J,Q)
            if flag:status[k,0]=4;break
            if it>=warm:
                sums[(it-warm)//(count//20)]+=logs;status[k,1]+=1
                soft_stats[k,0]+=peak;soft_stats[k,1]+=entropy;soft_stats[k,2]+=peak<.99
            if it>=warm+count-100000:states[k,it-(warm+count-100000)]=x
        if status[k,0]==0:rates[k]=sums/(count//20);soft_stats[k]/=count
    return rates,states,status,minima,soft_stats


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    with np.load(args.source) as f:
        ids=np.flatnonzero(f['flags']);params=f['params'][ids];points=f['points'][ids];counts=f['counts'][ids];kinds=f['kinds'][ids];base_ids=f['base_ids'][ids];relative=f['relative_tau'][ids]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),warm=30000,observed=500000,windows=20,initial=[.317,-.211,.419],selection='Every short-scout flag, including six hard-projection controls.',scope='Finite-time numerical experiment. Positive temperature replaces nearest projection by Gaussian posterior averaging; no novelty claim.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'search_soft_constellations.py','search_raar_constellations.py','follow_spherical_constellations.py','search_spherical_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[];begin=time.time()
    for offset in range(0,len(ids),8):
        sel=slice(offset,offset+8);rates,states,status,minima,stats=follow(params[sel],points[sel],counts[sel],30000,500000)
        for j,ident in enumerate(ids[sel]):
            k=offset+j
            np.savez_compressed(args.out/f'case-{ident}.npz',params=params[k],points=points[k,:counts[k]],kind=kinds[k],base_id=base_ids[k],relative_tau=relative[k],states=states[j],rates=rates[j],status=status[j],minima=minima[j],soft_stats=stats[j])
            row=dict(id=int(ident),base_id=int(base_ids[k]),relative_tau=float(relative[k]),beta=float(params[k,0]),tau=float(params[k,1]),status=status[j].tolist(),spectrum_qr_order=rates[j].mean(0).tolist(),minimum_leading_window=float(rates[j,:,0].min()),minimum_radius=float(minima[j,0]),minimum_peak_weight=float(minima[j,1]),soft_stats=stats[j].tolist())
            if status[j,0]==0:
                z=states[j];ev=np.linalg.eigvalsh(np.cov(z.T));row['variance_fractions']=(ev/ev.sum()).tolist();scale=max(float(np.linalg.norm(np.ptp(z,axis=0))),1.)
                residuals=[float(np.max(np.linalg.norm(z[-2048+k:]-z[-2048:-k],axis=1))) for k in range(1,513)];periods=[k+1 for k,v in enumerate(residuals) if v<1e-8*scale]
                row['detected_period']=min(periods) if periods else None;row['return_residuals']=residuals
            rows.append(row)
        (args.out/'partial.json').write_text(json.dumps(dict(protocol,cases=rows,seconds=time.time()-begin),indent=2)+'\n',encoding='utf-8',newline='\n')
        print(json.dumps(dict(completed=len(rows),total=len(ids),seconds=time.time()-begin)),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(protocol,cases=rows,seconds=time.time()-begin),indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

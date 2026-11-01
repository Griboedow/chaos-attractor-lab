"""Full local derivative spectra and fresh-start clouds for all scout flags."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_raar_constellations import step


@njit(cache=True,parallel=True)
def follow(params,points,counts,warm,count):
    n=len(params);rates=np.full((n,20,3),np.nan);states=np.full((n,100000,3),np.nan);status=np.zeros((n,2),np.int64);minima=np.full((n,2),np.inf)
    for k in prange(n):
        beta=params[k,0];p=points[k,:counts[k]];x=np.array([.317,-.211,.419]);Q=np.eye(2);sums=np.zeros((20,3))
        for it in range(warm+count):
            x,J,gap,axis,flag=step(x,beta,p);minima[k,0]=min(minima[k,0],gap);minima[k,1]=min(minima[k,1],axis)
            if flag or not np.isfinite(x).all():status[k,0]=flag if flag else 3;break
            R=J[:2,:2]@Q;a=np.sqrt(R[0,0]**2+R[1,0]**2)
            if a<1e-280:status[k,0]=4;break
            Q[0,0]=R[0,0]/a;Q[1,0]=R[1,0]/a;dot=Q[0,0]*R[0,1]+Q[1,0]*R[1,1]
            R[0,1]-=dot*Q[0,0];R[1,1]-=dot*Q[1,0];b=np.sqrt(R[0,1]**2+R[1,1]**2)
            if b<1e-280:status[k,0]=4;break
            Q[0,1]=R[0,1]/b;Q[1,1]=R[1,1]/b
            if it>=warm:
                window=(it-warm)//(count//20);sums[window,0]+=np.log(a);sums[window,1]+=np.log(b);sums[window,2]+=np.log(beta);status[k,1]+=1
            if it>=warm+count-100000:states[k,it-(warm+count-100000)]=x
        if status[k,0]==0:rates[k]=sums/(count//20)
    return rates,states,status,minima


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6)
    with np.load(args.source) as f:
        ids=np.flatnonzero(f['flags']);params=f['params'][ids];points=f['points'][ids];counts=f['counts'][ids];kinds=f['kinds'][ids]
    begin=time.time();rows=[]
    for offset in range(0,len(ids),6):
        sel=slice(offset,offset+6);rates,states,status,minima=follow(params[sel],points[sel],counts[sel],30000,500000)
        for j,i in enumerate(ids[sel]):
            k=offset+j;np.savez_compressed(args.out/f'case-{i}.npz',states=states[j],params=params[k],points=points[k,:counts[k]],rates=rates[j],status=status[j],minima=minima[j],kind=kinds[k])
            row=dict(id=int(i),kind=int(kinds[k]),status=status[j].tolist(),spectrum_xy_then_z=rates[j].mean(0).tolist(),minimum_leading_xy_window=float(rates[j,:,0].min()),minimum_gap=float(minima[j,0]),minimum_axis=float(minima[j,1]))
            if status[j,0]==0:
                ev=np.linalg.eigvalsh(np.cov(states[j].T));row['variance_fractions']=(ev/ev.sum()).tolist();scale=max(float(np.linalg.norm(np.ptp(states[j],axis=0))),1.)
                residuals=[float(np.max(np.linalg.norm(states[j,-2048+lag:]-states[j,-2048:-lag],axis=1))) for lag in range(1,513)]
                periods=[lag+1 for lag,v in enumerate(residuals) if v<1e-8*scale];row['detected_period']=min(periods) if periods else None;row['return_residuals']=residuals
            rows.append(row)
        print(json.dumps(dict(completed=min(offset+6,len(ids)),total=len(ids),seconds=time.time()-begin)),flush=True)
    report=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),warm=30000,observed=500000,initial=[.317,-.211,.419],cases=rows,scope='Local derivative spectrum of a discontinuous nearest-point map; no asymptotic or historical novelty proof. All original flags followed.',seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    for name in [Path(__file__).name,'search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())

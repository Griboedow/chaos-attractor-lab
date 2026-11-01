"""Fresh-start full spectra and native states for every flagged cyclic-RAAR scout."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_cyclic_raar import geometry,step


@njit(cache=True,parallel=True)
def follow(parameters,warm,count):
    rates=np.zeros((len(parameters),20,3));tails=np.empty((len(parameters),200000,3));axes=np.full(len(parameters),np.inf)
    for k in prange(len(parameters)):
        p=parameters[k];c,n,r=geometry(p);s=np.array([.317,-.211,.419]);Q=np.eye(3)
        for t in range(warm+count):
            s,j,axis,flag=step(s,p[0],c,n,r);axes[k]=min(axes[k],axis)
            if flag or not np.isfinite(s).all():raise ValueError('Invalid state')
            Q=j@Q
            for col in range(3):
                for prior in range(col):Q[:,col]-=np.dot(Q[:,prior],Q[:,col])*Q[:,prior]
                length=np.linalg.norm(Q[:,col])
                if length<1e-280:raise ValueError('Degenerate tangent')
                Q[:,col]/=length
                if t>=warm:rates[k,(t-warm)//(count//20),col]+=np.log(length)
            if t>=warm+count-200000:tails[k,t-(warm+count-200000)]=s
    return rates/(count//20),tails,axes


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6)
    with np.load(args.source) as f:
        ids=np.flatnonzero(f['flags']);p=f['params'][ids]
    begin=time.time();rates,states,axes=follow(p,args.warm,args.count)
    report=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),warm=args.warm,observed=args.count,initial=[.317,-.211,.419],units='per complete three-pair cycle',cases=[])
    for j,i in enumerate(ids):
        np.savez_compressed(args.out/f'case-{i}.npz',states=states[j],params=p[j],rates=rates[j],minimum_axis=axes[j])
        residuals=[float(np.max(np.linalg.norm(states[j,-2048+lag:]-states[j,-2048:-lag],axis=1))) for lag in range(1,513)]
        scale=max(float(np.linalg.norm(np.ptp(states[j],axis=0))),1.)
        periods=[lag+1 for lag,v in enumerate(residuals) if v<1e-8*scale]
        report['cases'].append(dict(id=int(i),spectrum=rates[j].mean(0).tolist(),minimum_leading_window=float(rates[j,:,0].min()),minimum_axis=float(axes[j]),detected_period=min(periods) if periods else None,return_residuals=residuals))
    report['seconds']=time.time()-begin
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    for name in [Path(__file__).name,'search_cyclic_raar.py','search_raar_circles.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps({**report,'cases':[{k:v for k,v in row.items() if k!='return_residuals'} for row in report['cases']]}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--warm',type=int,default=50000);p.add_argument('--count',type=int,default=1000000);main(p.parse_args())

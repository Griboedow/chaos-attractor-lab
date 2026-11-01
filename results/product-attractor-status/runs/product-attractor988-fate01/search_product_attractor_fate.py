"""Search for capture, persistent chaos and near-singular events at fixed Case 988.

The declared cohorts are not a basin-volume sample. Finite survival, including
positive tangent growth, is not an attracting-set certificate.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.stats import qmc


@njit(cache=True)
def fast_step(x,w,beta,tau,levels):
    m=np.empty(3);d=np.empty(3)
    for j in range(3):
        a=(levels[j,1]-levels[j,0])/2;c=(levels[j,1]+levels[j,0])/2
        h=2*a*(x[j]-c)/tau;t=np.tanh(h);m[j]=c+a*t
        e=np.exp(-2*abs(h));d[j]=8*a*a/tau*e/(1+e)**2
    v=2*m-x;r=np.sqrt(np.dot(v,v))
    if r<1e-14:return x,w,r,0.,1
    u=v/r;y=beta*x+(1-2*beta)*m+beta*u
    t=(2*d-1)*w;z=beta*w+(1-2*beta)*d*w+beta/r*(t-u*np.dot(u,t))
    length=np.sqrt(np.dot(z,z))
    if length<1e-280 or not np.isfinite(length):return y,w,r,0.,2
    return y,z/length,r,np.log(length),0


@njit(cache=True,parallel=True)
def follow(starts,beta,tau,levels,warm,count,block):
    n=len(starts);nb=count//block;tail=np.empty((n,2048,3));rates=np.zeros((n,nb))
    minima=np.full((n,nb),np.inf);minimum_states=np.zeros((n,3));minimum_times=np.zeros(n,np.int64)
    ends=np.empty((n,nb,3));flags=np.zeros(n,np.int64);bounds=np.zeros((n,2,3));bounds[:,0]=np.inf;bounds[:,1]=-np.inf
    for k in prange(n):
        x=starts[k].copy();w=np.array([1.,.371,-.217]);w/=np.linalg.norm(w);lowest=np.inf
        for it in range(warm+count):
            y,w,r,g,flag=fast_step(x,w,beta,tau,levels)
            if flag or not np.isfinite(y).all():flags[k]=flag or 3;break
            if it>=warm:
                t=it-warm;ib=t//block;rates[k,ib]+=g;minima[k,ib]=min(minima[k,ib],r)
                if r<lowest:lowest=r;minimum_states[k]=x;minimum_times[k]=t
                if (t+1)%block==0:ends[k,ib]=y
                for j in range(3):bounds[k,0,j]=min(bounds[k,0,j],y[j]);bounds[k,1,j]=max(bounds[k,1,j],y[j])
                if t>=count-2048:tail[k,t-(count-2048)]=y
            x=y
    return tail,rates/block,minima,minimum_states,minimum_times,ends,flags,bounds


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time();set_num_threads(4)
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    with np.load(args.source) as data:reference=data['states'];beta,tau=data['params'];levels=data['levels']
    cycle=np.load(args.cycle)['orbit'][0]
    z=qmc.Sobol(4,scramble=True,seed=2026110201).random_base2(6);v=2*z[:,:3]-1
    indices=np.floor(z[:,3]*len(reference)).astype(int)
    starts=np.vstack([2.2*v,reference[indices]+1e-5*v,reference[indices]+.01*v,cycle+.01*v])
    labels=['global_cube','reference_1e-5','reference_0.01','cycle_0.01']
    protocol=dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),cycle_sha256=hashlib.sha256(args.cycle.read_bytes()).hexdigest(),
        seed=2026110201,labels=labels,starts_per_group=64,starts=starts.tolist(),reference_indices=indices.tolist(),warm=args.warm,count=args.count,block=args.block,
        formula='Float64 binary tanh posterior; derivative uses a stable exponential sech-squared expression.',scope=__doc__)
    save=lambda name,v:(args.out/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    save('protocol.json',protocol);parts=[]
    for offset in range(0,len(starts),16):
        result=follow(starts[offset:offset+16],beta,tau,levels,args.warm,args.count,args.block)
        np.savez_compressed(args.out/f'block-{offset:03}.npz',starts=starts[offset:offset+16],tail=result[0],rates=result[1],minima=result[2],minimum_states=result[3],minimum_times=result[4],ends=result[5],flags=result[6],bounds=result[7])
        parts.append(result)
        progress=dict(completed=min(offset+16,len(starts)),total=len(starts),minimum_radius=float(result[2].min()),minimum_leading_block=float(result[1].min()),flags=int(np.count_nonzero(result[6])),seconds=time.time()-begin)
        save('progress.json',progress);print(json.dumps(progress),flush=True)
    arrays=[np.concatenate([v[k] for v in parts]) for k in range(8)]
    tail,rates,minima,near,when,ends,flags,bounds=arrays;periods=np.zeros(len(starts),int)
    for k in range(len(starts)):
        if flags[k]:continue
        for lag in list(range(1,129))+[256,512]:
            if np.max(np.linalg.norm(tail[k,-512:]-tail[k,-512-lag:-lag],axis=1))<1e-9:periods[k]=lag;break
    np.savez_compressed(args.out/'fate.npz',starts=starts,tail=tail,rates=rates,minima=minima,minimum_states=near,minimum_times=when,ends=ends,flags=flags,bounds=bounds,periods=periods,params=[beta,tau],levels=levels)
    groups=[]
    for i,label in enumerate(labels):
        sl=slice(i*64,(i+1)*64)
        groups.append(dict(label=label,starts=64,flags=int(np.count_nonzero(flags[sl])),all_blocks_positive=int(np.sum(rates[sl].min(1)>0)),detected_periods=periods[sl].tolist(),
            mean_leading_range=[float(rates[sl].mean(1).min()),float(rates[sl].mean(1).max())],minimum_radius=float(minima[sl].min())))
    save('report.json',dict(protocol,groups=groups,total_observed=len(starts)*args.count,seconds=time.time()-begin));print(json.dumps(groups),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','cycle','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--warm',type=int,default=20000);p.add_argument('--count',type=int,default=4000000);p.add_argument('--block',type=int,default=200000)
    main(p.parse_args())

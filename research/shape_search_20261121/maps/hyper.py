"""Search for structured three-state clouds with two positive finite-time rates."""
import argparse
import json
from pathlib import Path
from search import CACHE, START, design, shape_features
import numpy as np
from numba import njit, prange


@njit(cache=True,parallel=True)
def two_rates(params,burn,observed):
    rates=np.zeros((len(params),4,2));tails=np.empty((len(params),640,3))
    width=observed//4
    for h in prange(len(params)):
        p=params[h];z=START.copy()
        q=np.eye(3)[:2].copy()
        for n in range(burn+observed):
            for i in range(3):
                z[i]*=p[0]
                q[0,i]*=p[0];q[1,i]*=p[0]
            z[0]+=1-p[0]
            for i in range(3):
                j,k=(i+1)%3,(i+2)%3
                w=z[k];cub=p[7+i]*w*w;angle=p[1+i]+(p[4+i]+cub)*w
                co,si=np.cos(angle),np.sin(angle)
                u,y=co*z[i]-si*z[j],si*z[i]+co*z[j]
                for t in range(2):
                    da=(p[4+i]+3*cub)*q[t,k]
                    vi,vj=co*q[t,i]-si*q[t,j]-y*da,si*q[t,i]+co*q[t,j]+u*da
                    q[t,i],q[t,j]=vi,vj
                z[i],z[j]=u,y
            for j in range(2):
                if j==1:
                    dot=q[1,0]*q[0,0]+q[1,1]*q[0,1]+q[1,2]*q[0,2]
                    for i in range(3):
                        q[1,i]-=dot*q[0,i]
                norm=np.sqrt(q[j,0]*q[j,0]+q[j,1]*q[j,1]+q[j,2]*q[j,2])
                for i in range(3):
                    q[j,i]/=norm
                if n>=burn:
                    rates[h,(n-burn)//width,j]+=np.log(norm)/width
            if n>=burn+observed-640:
                tails[h,n-burn-observed+640]=z
    return rates,tails


def run(cache,batch,count):
    q=design(batch+800,count)
    q[:,0]=.66+.305*(q[:,0]-.34)/.54
    q[:,4:7]*=.90
    q[:,7:10]*=.80
    features=np.full((count,18),np.nan);rates=np.empty((count,4,2))
    for begin in range(0,count,128):
        pars=q[begin:begin+128]
        w,tails=two_rates(pars,3072,8192)
        rates[begin:begin+len(pars)]=w
        for j in np.flatnonzero((w[:,:,1].min(axis=1)>.009)&(w[:,:,0].mean(axis=1)<1.05)):
            try:
                features[begin+j]=shape_features(tails[j])
            except ValueError:
                pass
    out=cache/f"batch-{batch:03}";out.mkdir(parents=True,exist_ok=True)
    ids=10_000_000+100_000*batch+np.arange(count)
    np.savez_compressed(out/"scout.npz",case=ids,params=q,rates=rates[:,:,0],second_rates=rates[:,:,1],features=features)
    print(batch,count,int((rates[:,:,1].min(axis=1)>.009).sum()),int(np.isfinite(features[:,0]).sum()),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE/"hyper")
    p.add_argument("--first",type=int,default=0)
    p.add_argument("--batches",type=int,default=8)
    p.add_argument("--count",type=int,default=4096)
    a=p.parse_args()
    for b in range(a.first,a.first+a.batches):
        run(a.cache,b,a.count)

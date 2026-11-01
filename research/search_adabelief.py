"""Frozen 4101-lane exploratory design, including five declared controls.

The scalar arithmetic intentionally preserves the original exploratory driver.
Screens indicate finite-time growth only; use full-spectrum follow-up controls.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import qmc


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    u=qmc.Sobol(2,scramble=True,seed=20261101).random_base2(12)
    b=np.r_[.05+.94*u[:,0],.9,.9,.9,.5,.8]
    c=np.r_[.05+.9495*u[:,1],.999,.99,.9,.9,.99]
    n=len(b);x=np.full(n,.73);m=np.zeros(n);s=np.zeros(n)
    v=np.tile(np.array([1.,.3,.2])[:,None],(1,n));v/=np.linalg.norm(v,axis=0)
    windows=[];acc=np.zeros(n);lo=np.full(n,np.inf);hi=-lo;t0=time.time()
    for t in range(30000):
        mn=b*m+(1-b)*x;r=b*(x-m)
        sn=c*s+(1-c)*r*r+1e-12;root=np.sqrt(sn);q=root+1e-14
        ds=2*(1-c)*b*r*(v[0]-v[1])+c*v[2]
        dm=(1-b)*v[0]+b*v[1]
        dx=v[0]-dm/q+mn*ds/(2*root*q*q)
        v=np.array([dx,dm,ds]);norm=np.linalg.norm(v,axis=0);v/=norm
        x=x-mn/q;m=mn;s=sn
        if not np.isfinite(x+m+s).all():raise FloatingPointError(f'invalid lane at {t}')
        if t>=10000:
            acc+=np.log(norm);lo=np.minimum(lo,x);hi=np.maximum(hi,x)
            if (t+1)%5000==0:windows.append(acc/5000);acc[:]=0
    w=np.array(windows)
    np.savez_compressed(args.output/'pilot.npz',b=b,c=c,windows=w,end=np.array([x,m,s]).T,lo=lo,hi=hi)
    report=dict(seed=20261101,sobol_lanes=4096,additional_controls=5,warmup=10000,
                observation=20000,window=5000,initial=[.73,0,0],threshold=.005,
                positive_all_windows=int((w.min(axis=0)>.005).sum()),seconds=time.time()-t0)
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()

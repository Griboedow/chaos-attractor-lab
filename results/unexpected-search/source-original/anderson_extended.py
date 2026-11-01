"""Regularized scalar AA(2), including monotone bounded base maps."""
from pathlib import Path
import json,time
import numpy as np
from scipy.stats import qmc
from numba import njit,prange,set_num_threads
set_num_threads(4)

@njit
def g_and_derivative(x,k,c,kind):
    if kind==4:return c+k*x,k
    if kind==0:return c+k*np.sin(x),k*np.cos(x)
    if kind==1:return c+k*np.cos(x),-k*np.sin(x)
    if kind==2:
        a=np.tanh(x);return c+k*a,k*(1-a*a)
    return c+k*x/(1+x*x),k*(1-x*x)/(1+x*x)**2

@njit
def step(z,v,p):
    k,c,reg,kind=p
    g=np.empty(3);gp=np.empty(3)
    for j in range(3):g[j],gp[j]=g_and_derivative(z[j],k,c,int(kind))
    r=g-z;mu=r.mean();dev=r-mu;den=reg+np.dot(dev,dev)
    gm=g.mean();cross=np.dot(dev,g);nxt=gm-mu*cross/den
    grad=np.empty(3)
    for j in range(3):
        rp=gp[j]-1
        grad[j]=gp[j]/3-(rp*cross/3+mu*((g[j]-gm)*rp+dev[j]*gp[j]))/den+mu*cross*2*dev[j]*rp/(den*den)
    return np.array([z[1],z[2],nxt]),np.array([v[1],v[2],np.dot(grad,v)])

@njit(parallel=True)
def survey(params,starts,warm,observed,keep):
    count=len(params);windows=np.full((count,5),np.nan);tails=np.zeros((count,keep,3));ends=np.zeros((count,3))
    status=np.zeros(count,np.int64);steps=np.zeros(count,np.int64);bounds=np.zeros(count)
    for lane in prange(count):
        p=params[lane];z=starts[lane].copy();v=np.array([1.,0.,0.]);acc=0.
        for t in range(warm+observed):
            z,v=step(z,v,p);norm=np.linalg.norm(v)
            if not np.isfinite(z).all() or np.max(np.abs(z))>1e6 or not np.isfinite(norm) or norm==0:
                status[lane]=-1;break
            bounds[lane]=max(bounds[lane],np.max(np.abs(z)));v/=norm
            if np.max(z)-np.min(z)<1e-10 and abs(g_and_derivative(z[2],p[0],p[1],int(p[3]))[0]-z[2])<1e-10:
                status[lane]=0;break
            if t>=warm:
                acc+=np.log(norm)
                if (t-warm+1)%(observed//5)==0:
                    windows[lane,(t-warm+1)//(observed//5)-1]=acc/(observed//5);acc=0.
                offset=t-(warm+observed-keep)
                if offset>=0:tails[lane,offset]=z
        else:status[lane]=1
        steps[lane]=t+1;ends[lane]=z
    return windows,tails,ends,status,steps,bounds

def main():
    out=Path(__file__).parent/'anderson-v2';out.mkdir(exist_ok=False)
    u=qmc.Sobol(6,scramble=True,seed=2026110133).random_base2(12)
    p=np.empty((4096,4));p[:,1]=0.;p[:,2]=10**(-8+10*u[:,2])
    p[:1024,0]=10**(-.05+1.4*u[:1024,0]);p[:1024,1]=6*u[:1024,1]-3;p[:1024,3]=2
    p[1024:2048,0]=10**(2*u[1024:2048,0]);p[1024:2048,3]=3
    p[2048:3072,0]=(2*u[2048:3072,0]-1)*20;p[2048:3072,3]=4
    p[3072:,0]=.5+3.5*u[3072:,0];p[3072:,1]=2*u[3072:,1]-1;p[3072:,3]=0
    starts=(2*u[:,3:6]-1)*np.maximum(1,np.sqrt(p[:,2,None]))*2
    for i in range(4096):
        if i%1024<512:
            starts[i,1]=g_and_derivative(starts[i,0],p[i,0],p[i,1],int(p[i,3]))[0]
            starts[i,2]=g_and_derivative(starts[i,1],p[i,0],p[i,1],int(p[i,3]))[0]
    begin=time.time();w,z,end,status,steps,bounds=survey(p,starts,5000,20000,2000)
    ids=np.where((status==1)&(w.min(1)>.005))[0]
    np.savez_compressed(out/'data.npz',parameters=p,starts=starts,windows=w,end=end,status=status,steps=steps,max_abs=bounds,positive_ids=ids,positive_tails=z[ids])
    print(json.dumps({'seconds':time.time()-begin,'statuses':{str(int(i)):int((status==i).sum()) for i in np.unique(status)},'positives':ids.tolist(),'max_lle':float(np.nanmax(w.mean(1))) if np.isfinite(w).any() else None}),flush=True)

if __name__=='__main__':main()

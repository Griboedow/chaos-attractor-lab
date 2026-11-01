"""Exploratory fixed-data one-unit FastICA dynamics; no algorithm modification."""
from pathlib import Path
import argparse, hashlib, json, math, platform, time
import numpy as np
from numba import njit, prange

@njit(cache=True)
def step(w,v,X,kind,alpha):
    n,d=X.shape;h=np.zeros(d);dv=np.zeros(d);avgp=0.;avgppv=0.
    for i in range(n):
        u=0.;z=0.
        for j in range(d):u+=X[i,j]*w[j];z+=X[i,j]*v[j]
        if kind==0:
            g=math.tanh(alpha*u);gp=alpha*(1-g*g);gpp=-2*alpha*gp*g
        elif kind==1:
            e=math.exp(-u*u/2);g=u*e;gp=(1-u*u)*e;gpp=(u*u*u-3*u)*e
        else:g=u*u*u;gp=3*u*u;gpp=6*u
        avgp+=gp/n;avgppv+=gpp*z/n
        for j in range(d):h[j]+=X[i,j]*g/n;dv[j]+=X[i,j]*gp*z/n
    for j in range(d):h[j]-=avgp*w[j];dv[j]-=avgp*v[j]+avgppv*w[j]
    norm=np.sqrt(np.dot(h,h))
    if norm<1e-14:return w,v,norm,0.
    f=h/norm;dv-=f*np.dot(f,dv);dv/=norm
    gain=np.sqrt(np.dot(dv,dv))
    if gain<1e-200:return f,v,norm,1e-200
    return f,dv/gain,norm,gain

@njit(cache=True)
def run(X,w,kind,alpha,warmup,steps,window):
    w=w/np.sqrt(np.dot(w,w));v=np.ones(len(w));v-=w*np.dot(w,v);v/=np.sqrt(np.dot(v,v))
    rates=np.zeros(steps//window);minimum=1e9;tail=np.empty((512,len(w)));status=0
    for i in range(warmup+steps):
        w,v,norm,gain=step(w,v,X,kind,alpha)
        minimum=min(minimum,norm)
        if norm<1e-14:status=1;break
        if i>=warmup:rates[(i-warmup)//window]+=math.log(gain)/window
        if i>=warmup+steps-512:tail[i-(warmup+steps-512)]=w
    period=0;best=1e9
    if status==0:
        for lag in range(1,129):
            error=0.
            for i in range(256,512):
                da=0.;db=0.
                for j in range(len(w)):da+=(tail[i,j]-tail[i-lag,j])**2;db+=(tail[i,j]+tail[i-lag,j])**2
                error=max(error,math.sqrt(min(da,db)))
            best=min(best,error)
            if error<1e-8 and period==0:period=lag
    return rates,minimum,period,best,status,w

@njit(parallel=True,cache=True)
def batch(X,counts,starts,kinds,alphas,warmup,steps,window):
    result=np.empty((len(X),steps//window+4));ends=np.empty_like(starts)
    for k in prange(len(X)):
        rates,mn,period,error,status,w=run(X[k,:counts[k]],starts[k],kinds[k],alphas[k],warmup,steps,window)
        result[k,:len(rates)]=rates;result[k,-4:]=np.array([mn,period,error,status]);ends[k]=w
    return result,ends

def dataset(seed,n,d,family):
    rng=np.random.default_rng(seed)
    if family==0:a=rng.normal(size=(n,d))
    elif family==1:
        a=rng.normal(size=(n,d));a/=np.linalg.norm(a,axis=1)[:,None]
    elif family==2:
        centers=rng.normal(size=(max(2,n//4),d))*2
        a=centers[rng.integers(len(centers),size=n)]+rng.normal(size=(n,d))*.3
    elif family==3:a=rng.standard_t(2.5,size=(n,d))
    elif family==4:a=rng.uniform(-1,1,size=(n,d))
    else:
        a=rng.normal(size=(n//2,d));a=np.concatenate([a,-a])
    a-=a.mean(axis=0);e,U=np.linalg.eigh(a.T@a/n)
    X=a@((U/np.sqrt(e))@U.T)
    return X,rng.normal(size=d)

def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False)
    protocol={'seed':2026110101,'cases':args.cases,'dimension':4,'sample_counts':[8,12,16,24,32,64],'families':['normal','sphere','Gaussian mixture','Student t(2.5)','uniform cube','antipodal normal'],'nonlinearities':['tanh(alpha*u), alpha in [1,2]','u*exp(-u*u/2)','u**3'],'warmup':1000,'observed':4000,'window':1000,'period_test':'Tail 512, lags 1..128, sup norm Euclidean up to sign, threshold 1e-8.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':platform.python_version(),'numpy':np.__version__}
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    rng=np.random.default_rng(protocol['seed']);X=np.zeros((args.cases,64,4));starts=np.empty((args.cases,4));ns=np.empty(args.cases,dtype=np.int64);kind=np.empty_like(ns);alpha=np.empty(args.cases);cases=[]
    for k in range(args.cases):
        n=int(rng.choice(protocol['sample_counts']));family=k%6;seed=int(rng.integers(1,2**32));ns[k]=n;kind[k]=(k//6)%3;alpha[k]=rng.uniform(1,2)
        x,w=dataset(seed,n,4,family);X[k,:n]=x;starts[k]=w;cases.append({'case':k,'seed':seed,'n':n,'family':family,'kind':int(kind[k]),'alpha':float(alpha[k]),'input_sha256':hashlib.sha256(x.tobytes()).hexdigest(),'mean_error':float(abs(x.mean(axis=0)).max()),'covariance_error':float(abs(x.T@x/n-np.eye(4)).max())})
    np.savez_compressed(out/'inputs.npz',X=X,counts=ns,starts=starts,kinds=kind,alphas=alpha)
    (out/'cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    start=time.time();stats,ends=batch(X,ns,starts,kind,alpha,1000,4000,1000)
    np.savez_compressed(out/'scout.npz',stats=stats,endpoints=ends)
    rows=[]
    for case,st in zip(cases,stats):rows.append({**case,'rates':st[:4].tolist(),'minimum_update_norm':st[-4],'period_up_to_sign':int(st[-3]),'best_recurrence_error':st[-2],'status':int(st[-1])})
    report={'wall_seconds':time.time()-start,'cases':rows,'retained':[v['case'] for v in rows if min(v['rates'])>.005 and v['period_up_to_sign']==0 and v['minimum_update_norm']>1e-8 and v['status']==0]}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'wall_seconds':report['wall_seconds'],'retained':report['retained'],'top':sorted(rows,key=lambda v:min(v['rates']),reverse=True)[:10]}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--cases',type=int,default=576);main(p.parse_args())

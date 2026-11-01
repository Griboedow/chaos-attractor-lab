"""Reproducible long-run and independent controls for scalar AdaBelief.

Run from any directory. Output must not already exist. Optional Numba compilation
uses no fast-math, no noise, no clipping, and no state/tangent resets in a run.
"""
import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np
from numba import njit
from adabelief_map import step, jacobian, inverse_step, spectrum


@njit(cache=False)
def compiled_run(initial, b, c, eps_s, eps_out, warmup, steps, window, keep,
                 corrected_warmup=False):
    z=initial.copy(); q=np.eye(3)
    windows=np.zeros((steps//window,3)); acc=np.zeros(3)
    tail=np.zeros((keep,3)); low=np.full(3,np.inf); high=-low
    moments=np.zeros((steps//window,4)); moment=np.zeros(4)
    for t in range(warmup+steps):
        x,m,s=z
        mn=b*m+(1-b)*x
        residual=x-mn
        sn=c*s+(1-c)*residual*residual+eps_s
        root=np.sqrt(sn); denom=root+eps_out
        j=np.zeros((3,3))
        j[1,0]=1-b; j[1,1]=b
        j[2,0]=2*(1-c)*b*residual
        j[2,1]=-j[2,0]; j[2,2]=c
        for k in range(3):
            j[0,k]=-j[1,k]/denom+mn*j[2,k]/(2*root*denom*denom)
        j[0,0]+=1
        a=j@q
        growth=np.zeros(3)
        # Modified Gram-Schmidt with reorthogonalization, every single iterate.
        for k in range(3):
            for repeat in range(2):
                for h in range(k):
                    dot=0.
                    for i in range(3): dot+=a[i,k]*q[i,h]
                    for i in range(3): a[i,k]-=dot*q[i,h]
            norm=np.sqrt(np.sum(a[:,k]**2))
            growth[k]=np.log(norm)
            q[:,k]=a[:,k]/norm
        if corrected_warmup and t<warmup:
            z[0]=x-(mn/(1-b**(t+1)))/(np.sqrt(sn/(1-c**(t+1)))+eps_out)
        else:
            z[0]=x-mn/denom
        z[1]=mn;z[2]=sn
        if not np.isfinite(z).all() or sn<=0:
            raise ValueError('invalid state')
        if t>=warmup:
            acc+=growth
            low=np.minimum(low,z);high=np.maximum(high,z)
            moment+=np.array([z[0],z[0]**2,z[1]**2,z[2]])
            n=t-warmup+1
            if n%window==0:
                windows[n//window-1]=acc/window;acc[:]=0
                moments[n//window-1]=moment/window;moment[:]=0
            if n>steps-keep: tail[n-(steps-keep)-1]=z
    return z,windows,tail,low,high,moments


def run(initial, b, c, *, eps_s=1e-12, eps_out=1e-14, warmup=20000,
        steps=1000000, window=10000, keep=100000, corrected_warmup=False):
    if steps%window or keep>steps or keep<0: raise ValueError('invalid run lengths')
    result=compiled_run(np.array(initial,dtype=float),b,c,eps_s,eps_out,
                        warmup,steps,window,keep,corrected_warmup)
    names=['end','windows','tail','minimum','maximum','moments']
    data=dict(zip(names,result))
    data['initial']=np.array(initial);data['parameters']=np.array([b,c,eps_s,eps_out])
    data['lengths']=np.array([warmup,steps,window,keep])
    data['corrected_warmup']=corrected_warmup
    data['spectrum']=data['windows'].mean(axis=0)
    return data


def brief(data):
    b,c=data['parameters'][:2]
    return {key:data[key].tolist() for key in ['initial','parameters','lengths','spectrum','minimum','maximum']} | {
        'volume_error':float(abs(data['spectrum'].sum()-np.log(b*c))),
        'window_min':data['windows'].min(axis=0).tolist(),
        'window_max':data['windows'].max(axis=0).tolist(),
        'corrected_warmup':bool(data['corrected_warmup'])}


def controls():
    rng=np.random.default_rng(20261102)
    states=rng.uniform(-1,1,(100,3));states[:,2]=10**rng.uniform(-4,1,100)
    maximum_complex=0.;maximum_inverse=0.;maximum_det=0.
    for b,c in [(.52,.15),(.9,.99),(.9,.999),(.2,.6)]:
        analytic=jacobian(states,b,c)
        numerical=np.empty_like(analytic)
        for k in range(3):
            zz=states.astype(complex);zz[:,k]+=1e-24j
            numerical[:,:,k]=step(zz,b,c).imag/1e-24
        maximum_complex=max(maximum_complex,float(np.max(abs(analytic-numerical)/(1+abs(numerical)))))
        maximum_inverse=max(maximum_inverse,float(abs(inverse_step(step(states,b,c),b,c)-states).max()))
        maximum_det=max(maximum_det,float(abs(np.linalg.det(analytic)-b*c).max()))
    # Householder QR / NumPy trajectory versus independently assembled compiled MGS.
    reference=spectrum([[.73,0,0]],.52,.15,warmup=1000,steps=20000,window=5000,keep=500)
    candidate=run([.73,0,0],.52,.15,warmup=1000,steps=20000,window=5000,keep=500)
    # Short trajectory comparisons only: chaotic tails need not be pointwise equal.
    z=np.array([.73,0.,0.]); short=[]
    for t in range(25): z=step(z,.52,.15);short.append(z)
    fast=run([.73,0,0],.52,.15,warmup=0,steps=25,window=25,keep=25)
    discrepancy=float(abs(fast['tail']-short).max())
    # Explicit sign-alternating period-two orbit in the exact same map.
    b,c=.1,.5; d=1e-14;e=1e-12
    k=(1-b)/(1+b);a=np.sqrt(((k/2-d)**2-e/(1-c)))/(2*b/(1+b))
    cycle=np.array([a,-k*a,(k/2-d)**2])
    next_state=step(cycle,b,c)
    multiplier=np.linalg.eigvals(jacobian(next_state,b,c)@jacobian(cycle,b,c))
    cycle_residual=float(abs(step(next_state,b,c)-cycle).max())
    result=dict(complex_step_relative_error=maximum_complex,inverse_error=maximum_inverse,
                determinant_error=maximum_det,short_trajectory_error=discrepancy,
                numpy_spectrum=reference['spectrum'][0].tolist(),
                compiled_spectrum=candidate['spectrum'].tolist(),
                cycle_state=cycle.tolist(),cycle_residual=cycle_residual,
                cycle_multiplier_moduli=abs(multiplier).tolist())
    assert maximum_complex<1e-10 and maximum_inverse<1e-10 and maximum_det<1e-9
    assert discrepancy<1e-10 and cycle_residual<1e-10
    assert abs(candidate['spectrum'].sum()-np.log(.52*.15))<1e-10
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mode',choices=['controls','long','ensemble','neighborhood','epsilon'],required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    start=time.time();reports=[]
    if args.mode=='controls':
        report=controls()
    else:
        jobs=[]
        if args.mode=='long':
            for b,c in [(.52,.15),(.9,.99),(.9,.999),(.9,.9),(.1,.5)]:
                jobs.append(dict(initial=[.73,0,0],b=b,c=c,steps=2000000,keep=200000))
        elif args.mode=='ensemble':
            rng=np.random.default_rng(20261103)
            cold=np.c_[np.r_[-10.,-1.,-.1,-.001,.001,.1,1.,10.],np.zeros((8,2))]
            general=np.c_[rng.uniform(-2,2,16),rng.uniform(-1,1,16),10**rng.uniform(-5,1,16)]
            for b,c in [(.52,.15),(.9,.99)]:
                for z in np.r_[cold,general]:
                    jobs.append(dict(initial=z.tolist(),b=b,c=c,steps=200000,keep=2048))
                for z in cold:
                    jobs.append(dict(initial=z.tolist(),b=b,c=c,steps=200000,keep=2048,corrected_warmup=True))
        elif args.mode=='neighborhood':
            for b,c in [(.52,.15),(.9,.99)]:
                for db in [-.002,-.001,0,.001,.002]:
                    for dc in [-.002,-.001,0,.001,.002]:
                        jobs.append(dict(initial=[.73,0,0],b=b+db,c=c+dc,steps=100000,keep=2048))
        elif args.mode=='epsilon':
            for b,c in [(.52,.15),(.9,.99),(.9,.999)]:
                for epsilon in [1e-16,1e-14,1e-12,1e-10,1e-8]:
                    jobs.append(dict(initial=[.73,0,0],b=b,c=c,eps_s=epsilon/1e-4,
                                     eps_out=epsilon/.01,steps=200000,keep=2048,corrected_warmup=True))
        for i,job in enumerate(jobs):
            data=run(**job);np.savez_compressed(args.output/f'run-{i:03d}.npz',**data)
            row=brief(data)|dict(id=i);reports.append(row)
            print(json.dumps(dict(mode=args.mode,id=i,spectrum=row['spectrum'])),flush=True)
        report={'runs':reports}
    report|={'mode':args.mode,'elapsed_seconds':time.time()-start,'python':platform.python_version(),
             'numpy':np.__version__,'platform':platform.platform()}
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report if args.mode=='controls' else {'runs':len(reports),'seconds':time.time()-start},indent=2),flush=True)


if __name__=='__main__': main()

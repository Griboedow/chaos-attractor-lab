"""RK23 controller dynamics for q'=p, p'=-q-2*zeta*p, in original axes.

The ODE is a stable linear damped oscillator for zeta>0. This accepted-step
map uses SciPy's componentwise weighted RMS error and standard factor rules.
State=(q/atol,p/atol,log proposed h); parameters=(zeta,rtol).
The finite floating-point clock and final-time clipping are omitted here.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.stats import qmc
from scipy.integrate import RK23
import scipy

@njit(cache=True)
def matrices(zeta):
    A=np.array([[0.,1.],[-1.,-2*zeta]])
    A2=A@A;A3=A2@A;A4=A3@A
    return A,A2,A3,A4

@njit(cache=True)
def oscillator_step(state,rtol,A,A2,A3,A4):
    u=state[:2].copy();x=state[2]
    du=np.array([[1.,0.,0.],[0.,1.,0.]])
    dx=np.array([0.,0.,1.]);rejected=False
    for attempt in range(100):
        h=np.exp(x);h2=h*h;h3=h2*h;h4=h3*h
        R=np.eye(2)+h*A+h2*A2/2+h3*A3/6
        E=(h3*A3+h4*A4)/48
        dr=h*A+h2*A2+h3*A3/2
        de=(3*h3*A3+4*h4*A4)/48
        new=R@u;error=E@u
        dnew=R@du+np.outer(dr@u,dx)
        derror=E@du+np.outer(de@u,dx)
        scale=np.empty(2);dscale=np.empty((2,3))
        for i in range(2):
            if abs(new[i])>abs(u[i]):
                scale[i]=1+rtol*abs(new[i]);dscale[i]=rtol*np.sign(new[i])*dnew[i]
            else:
                scale[i]=1+rtol*abs(u[i]);dscale[i]=rtol*np.sign(u[i])*du[i]
        ratio=error/scale
        norm=np.sqrt(np.sum(ratio*ratio)/2)
        if norm==0:
            logfactor=np.log(10.);df=np.zeros(3)
        else:
            derivative=np.zeros(3)
            for i in range(2):derivative+=ratio[i]*(derror[i]/scale[i]-error[i]/(scale[i]*scale[i])*dscale[i])
            derivative/=np.sum(ratio*ratio)
            logfactor=np.log(.9)-np.log(norm)/3;df=-derivative/3
        if norm<1:
            ceiling=0. if rejected else np.log(10.)
            if logfactor>ceiling:logfactor=ceiling;df*=0
            result=np.array([new[0],new[1],x+logfactor])
            jac=np.empty((3,3));jac[:2]=dnew;jac[2]=dx+df
            return result,jac,attempt,0
        if logfactor<np.log(.2):logfactor=np.log(.2);df*=0
        x+=logfactor;dx+=df;rejected=True
    return state,np.eye(3),100,1

@njit(cache=True,parallel=True)
def screen(params,warm,steps):
    n=len(params);rates=np.full((n,10),np.nan);bounds=np.empty((n,2,3));ends=np.empty((n,3));counts=np.zeros((n,3),dtype=np.int64);tails=np.full((n,2048,3),np.nan)
    for k in prange(n):
        A,A2,A3,A4=matrices(params[k,0]);rtol=params[k,1];s=np.array([1.,0.,np.log(.1)])
        v=np.array([1.,.371,-.217]);v/=np.sqrt(np.sum(v*v));r=np.zeros(10);lo=np.full(3,np.inf);hi=np.full(3,-np.inf)
        for t in range(warm+steps):
            if not np.all(np.isfinite(s)) or max(abs(s[0]),abs(s[1]))>1e10 or abs(s[2])>50:
                counts[k,0]=2;break
            s,j,reject,status=oscillator_step(s,rtol,A,A2,A3,A4)
            if status:counts[k,0]=1;break
            v=j@v;norm=np.sqrt(np.sum(v*v))
            if not np.isfinite(norm) or norm<1e-290:counts[k,0]=3;break
            v/=norm
            if t>=warm:
                r[(t-warm)//(steps//10)]+=np.log(norm);lo=np.minimum(lo,s);hi=np.maximum(hi,s);counts[k,1]+=1;counts[k,2]+=reject
            if t>=warm+steps-2048:tails[k,t-(warm+steps-2048)]=s
        if counts[k,0]==0:rates[k]=r/(steps//10)
        ends[k]=s;bounds[k,0]=lo;bounds[k,1]=hi
    return rates,bounds,ends,counts,tails

def verify():
    rng=np.random.default_rng(2026110174);diff=[];upstream=[]
    for _ in range(256):
        zeta=rng.uniform(.01,.99);rtol=10**rng.uniform(-12,-1);ms=matrices(zeta)
        s=np.r_[rng.uniform(-2,2,2),rng.uniform(-1.5,1.5)]
        new,jac,_,status=oscillator_step(s,rtol,*ms);assert status==0
        fd=np.column_stack([(oscillator_step(s+np.eye(3)[j]*1e-6,rtol,*ms)[0]-oscillator_step(s-np.eye(3)[j]*1e-6,rtol,*ms)[0])/2e-6 for j in range(3)])
        diff.append(float(np.max(abs(jac-fd))))
        A=ms[0];solver=RK23(lambda t,y:A@y,0.,s[:2],1e10,rtol=rtol,atol=1.,first_step=np.exp(s[2]))
        solver.step();expected=np.r_[solver.y,np.log(solver.h_abs)]
        upstream.append(float(np.max(abs(new-expected))))
    assert max(diff)<1e-6,max(diff)
    assert max(upstream)<1e-12,max(upstream)
    return dict(checks=256,max_finite_difference_jacobian_error=max(diff),max_unmodified_scipy_one_step_residual=max(upstream))

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(12);checks=verify()
    sobol=qmc.Sobol(2,scramble=True,seed=2026110175).random_base2(args.power)
    params=np.column_stack([.001+.998*sobol[:,0],10**(-12+11*sobol[:,1])])
    params=np.vstack([params,[[.2,1e-3],[.5,1e-3],[.8,1e-3]]]);begin=time.time()
    rates,bounds,ends,counters,tails=screen(params,20000,50000)
    flags=(counters[:,0]==0)&(rates.min(1)>.001)
    np.savez_compressed(args.out/'screen.npz',params=params,rates=rates,bounds=bounds,ends=ends,counters=counters,tails=tails,flags=flags)
    report=dict(scipy_version=scipy.__version__,checks=checks,count=len(params),warm=20000,observed=50000,windows=10,
                retained=int(flags.sum()),completed=int(np.sum(counters[:,0]==0)),seconds=time.time()-begin,
                seed=2026110175,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Controller state for a stable linear oscillator. Numerical behavior near atol, not physical ODE chaos.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--power',type=int,default=10);main(p.parse_args())

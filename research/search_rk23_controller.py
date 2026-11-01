"""SciPy-style RK23 accepted-step amplitude map on stable complex linear ODEs.

The reduction removes the S1 phase symmetry; any added phase exponent is zero.
It transcribes the installed SciPy acceptance/factor rule, omitting finite-time
endpoint clipping and floating-point t+h subtraction. An upstream check is
required before describing a candidate as reproduced by SciPy itself.
"""
import argparse,hashlib,inspect,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.integrate._ivp.rk import RK23,rk_step,RungeKutta
from scipy.stats import qmc
import scipy

@njit(cache=True)
def rk23_step(state,params):
    # state=(log proposed h, log |y|/atol); params=(angle,rtol).
    x,y=state;angle,rtol=params
    dx=np.array([1.,0.]);dy=np.array([0.,1.])
    lam=complex(-np.cos(angle),np.sin(angle))
    rejected=False
    for attempt in range(100):
        h=np.exp(x);z=h*lam
        R=1+z+z*z/2+z*z*z/6
        E=z*z*z*(1+z)/48
        modulus=abs(R)
        if modulus==0:
            # Exact arrival at y=0 leaves this positive-amplitude chart.
            return state,np.eye(2),0.,attempt,2
        mu=(z*(1+z+z*z/2)/R).real
        scale_term=rtol*np.exp(y)*max(1.,modulus)
        dscale=dy.copy()
        if modulus>1:dscale+=mu*dx
        if E==0:
            # The capped-factor branch is constant in a neighborhood of this
            # zero. Do not form E'/E before selecting that branch.
            logerror=-np.inf;de=np.zeros(2)
        else:
            nu=(3+z/(1+z)).real
            logerror=np.log(abs(E))+y-np.log1p(scale_term)
            de=nu*dx+dy-scale_term/(1+scale_term)*dscale
        logfactor=np.log(.9)-logerror/3
        df=-de/3
        if logerror<0:
            ceiling=0. if rejected else np.log(10.)
            if logfactor>ceiling:logfactor=ceiling;df*=0
            result=np.array([x+logfactor,y+np.log(modulus)])
            jac=np.empty((2,2));jac[0]=dx+df;jac[1]=dy+mu*dx
            return result,jac,np.angle(R),attempt,0
        if logfactor<np.log(.2):logfactor=np.log(.2);df*=0
        x+=logfactor;dx+=df;rejected=True
    return state,np.eye(2),0.,100,1

@njit(cache=True,parallel=True)
def screen(params,warm,steps):
    n=len(params);rates=np.full((n,10),np.nan);bounds=np.empty((n,2,2));ends=np.empty((n,2));counts=np.zeros((n,3),dtype=np.int64);states=np.full((n,1024,2),np.nan)
    for k in prange(n):
        s=np.array([np.log(.1),0.]);v=np.array([.9,.435889894]);total=np.zeros(10);lo=np.full(2,np.inf);hi=np.full(2,-np.inf)
        for t in range(warm+steps):
            if not np.all(np.isfinite(s)) or max(abs(s[0]),abs(s[1]))>100:
                counts[k,0]=2;break
            s,j,phase,reject,status=rk23_step(s,params[k])
            if status:counts[k,0]=1;break
            v=j@v;norm=np.sqrt(np.sum(v*v))
            if not np.isfinite(norm) or norm<1e-290:counts[k,0]=3;break
            v/=norm
            if t>=warm:
                total[(t-warm)//(steps//10)]+=np.log(norm);lo=np.minimum(lo,s);hi=np.maximum(hi,s);counts[k,1]+=1;counts[k,2]+=reject
            if t>=warm+steps-1024:states[k,t-(warm+steps-1024)]=s
        if counts[k,0]==0:rates[k]=total/(steps//10)
        bounds[k,0]=lo;bounds[k,1]=hi;ends[k]=s
    return rates,bounds,ends,counts,states

def verify():
    rng=np.random.default_rng(2026110172);errors=[];polynomial_errors=[]
    for _ in range(512):
        angle=rng.uniform(0,np.pi/2-.001);rtol=10**rng.uniform(-12,-1);p=np.array([angle,rtol]);s=rng.uniform(-1.5,1.5,2)
        _,j,_,_,status=rk23_step(s,p);assert status==0
        finite=np.column_stack([(rk23_step(s+np.eye(2)[i]*1e-6,p)[0]-rk23_step(s-np.eye(2)[i]*1e-6,p)[0])/2e-6 for i in range(2)])
        errors.append(float(np.max(abs(j-finite))))
        lam=-np.cos(angle)+1j*np.sin(angle);h=np.exp(s[0]);z=h*lam;y=np.array([np.exp(s[1])+0j]);K=np.zeros((4,1),dtype=complex)
        new,_=rk_step(lambda t,u:lam*u,0.,y,lam*y,h,RK23.A,RK23.B,RK23.C,K)
        expected=(1+z+z*z/2+z*z*z/6)*y
        polynomial_errors.append(float(max(abs(new-expected).max(),abs(h*(K.T@RK23.E)-z**3*(1+z)/48*y).max())))
    assert max(errors)<1e-6,max(errors)
    assert max(polynomial_errors)<1e-12,max(polynomial_errors)
    value,jac,_,_,status=rk23_step(np.zeros(2),np.array([0.,1e-3]))
    assert status==0 and np.max(abs(value-np.log([10.,1/3])))<1e-14
    assert np.max(abs(jac-np.array([[1.,0.],[-1.5,1.]])))<1e-14
    return dict(finite_difference_count=512,max_jacobian_error=max(errors),max_installed_rk_stage_polynomial_error=max(polynomial_errors),exact_zero_error_branch_check=True)

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);checks=verify()
    s=qmc.Sobol(2,scramble=True,seed=2026110173).random_base2(args.power)
    p=np.column_stack([s[:,0]*(np.pi/2-1e-4),10**(-14+13*s[:,1])])
    p=np.vstack([p,[[0,1e-3],[np.pi/4,1e-3],[np.pi/2-.01,1e-3]]])
    begin=time.time();rates,bounds,ends,counts,states=screen(p,20000,50000)
    flags=(counts[:,0]==0)&(rates.min(1)>.001)
    np.savez_compressed(args.out/'screen.npz',params=p,rates=rates,bounds=bounds,ends=ends,counters=counts,tails=states,flags=flags)
    report=dict(scipy_version=scipy.__version__,checks=checks,count=len(p),warm=20000,observed=50000,windows=10,
                positive=int(flags.sum()),completed=int(np.sum(counts[:,0]==0)),seconds=time.time()-begin,
                seed=2026110173,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Amplitude/step quotient. Phase symmetry contributes a neutral exponent, not an independent chaotic mechanism.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    (args.out/'installed-step-implementation.py').write_text(inspect.getsource(RungeKutta._step_impl),encoding='utf-8')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--power',type=int,default=10);main(p.parse_args())

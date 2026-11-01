"""Declared Yogi shape search; no historical uniqueness claim."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from numba import set_num_threads
from yogi_map import scout,step_tangent


def controls():
    rng=np.random.default_rng(2026110140);errors=[];deterrors=[]
    for _ in range(400):
        p=np.r_[rng.uniform(.02,.99,2),10**rng.uniform(-6,-1)];z=np.r_[rng.uniform(-2,2,2),rng.uniform(.02,3)]
        if abs(z[2]-z[0]**2)<1e-4:continue
        sign=np.sign(z[2]-z[0]**2);b,c,e=p
        def scalar(w):
            x,m,v=w;mn=b*m+(1-b)*x;vn=v-(1-c)*sign*x*x
            return np.array([x-mn/(np.sqrt(vn)+e),mn,vn])
        val,j=step_tangent(z,np.eye(3),p)
        cs=np.column_stack([np.imag(scalar(z.astype(complex)+1e-25j*np.eye(3)[k]))/1e-25 for k in range(3)])
        error=float(np.max(abs(j-cs))/(1+np.max(abs(j))));assert error<1e-12
        assert np.max(abs(val-scalar(z)))<1e-12
        deterror=float(abs(np.linalg.det(j)-b));assert deterror<1e-10
        errors.append(error);deterrors.append(deterror)
    return dict(samples=len(errors),maximum_scaled_jacobian_error=max(errors),maximum_determinant_error=max(deterrors),
        method='Independent complex-step evaluation within each smooth sign branch. The branch determinant is beta1. Discontinuity boundaries are excluded.')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    seed=2026110141;u=qmc.Sobol(3,scramble=True,seed=seed).random_base2(12)
    base=np.c_[.02+.97*u[:,0],.02+.97*u[:,1],10**(-6+5*u[:,2])]
    base[2048:,1]=1-10**(-4+3.95*u[2048:,1]);base[::2,2]=.1
    base=np.r_[base,[[.9,.999,.1],[.9,.999,1e-6],[.5,.5,.1],[.5,.5,1e-6],[.1,.1,.1]]]
    p=np.repeat(base,2,axis=0);starts=np.tile([.73,0.,0.],(len(p),1));starts[1::2,2]=.73**2
    protocol=dict(seed=seed,count=len(p),alpha=.01,objective='theta^2/2',coordinates='x=theta/alpha,m=M/alpha,v=V/alpha^2',
        warmup=10000,observed=20000,window=5000,tail=512,
        design='4096 Sobol parameter rows plus five controls, each with two prescribed initial second moments. beta1=.02+.97*u; first 2048 beta2=.02+.97*v, remaining beta2=1-10^(-4+3.95*v). Even parameter rows have e=.1; odd rows have e=10^(-6+5*w).',
        initialization='x=.73,m=0; paired v=0 and v=x^2. Autonomous Algorithm2 throughout, no bias correction.',
        primary_rule='Completed; all four leading derivative windows>.005; final-tail x range>1e-4; no tail return at lags1..32 under tolerance1e-9*(1+maximum state magnitude).',
        source='Zaheer et al2018 Algorithm2; source PDF retrieved before implementation. Optax scale_by_yogi inspected separately; it uses bias correction and positive initial accumulator, so this scout does not claim exact cold-start library equality.',
        scope='Discontinuous map; smooth-branch derivative growth is a screening statistic, not an attractor or novelty proof.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');(out/'controls.json').write_text(json.dumps(controls(),indent=2)+'\n')
    for name in ['search_yogi.py','yogi_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    before=time.perf_counter();windows,tails,ends,status,extrema,margins,changes=scout(p,starts,10000,20000,5000,512)
    periods=np.zeros(len(p),dtype=int)
    for i in range(len(p)):
        if status[i]!=0:continue
        tol=1e-9*(1+np.max(abs(tails[i])))
        for lag in range(1,33):
            if np.max(abs(tails[i,lag:]-tails[i,:-lag]))<tol:periods[i]=lag;break
    file=out/'scout.npz';np.savez_compressed(file,parameters=p,starts=starts,windows=windows,tails=tails,ends=ends,status=status,extrema=extrema,minimum_margins=margins,branch_changes=changes,periods=periods)
    keep=(status==0)&(windows.min(1)>.005)&(np.ptp(tails[:,:,0],axis=1)>1e-4)&(periods==0)
    result=dict(completed=int((status==0).sum()),declared=len(p),primary_count=int(keep.sum()),primary_ids=np.flatnonzero(keep).tolist(),
        zero_tail_count=int(np.sum(np.max(abs(tails[:,:,:2]),axis=(1,2))<1e-100)),
        controls=[dict(id=int(i),parameters=p[i].tolist(),start=starts[i].tolist(),windows=windows[i].tolist(),period=int(periods[i])) for i in range(8192,len(p))],
        seconds=time.perf_counter()-before,sha256=hashlib.sha256(file.read_bytes()).hexdigest())
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['completed','declared','primary_count','zero_tail_count','controls','seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

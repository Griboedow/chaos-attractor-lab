"""Declared paired-start scout for the existing TensorFlow Addons tanh Yogi."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from numba import set_num_threads
from yogi_tanh_map import scout,step_tangent


def controls():
    rng=np.random.default_rng(2026110161);errors=[];dets=[];values=[]
    for i in range(512):
        b,c=rng.uniform(.02,.995,2);e=10**rng.uniform(-6,-1);k=10**rng.uniform(-3,5)
        x=rng.uniform(-1,1);m=rng.uniform(-1,1)
        v=max(x*x+rng.uniform(-2,2)/k,1e-5) if i%2 else 10**rng.uniform(-5,1)
        z=np.array([x,m,v]);p=np.array([b,c,e,k]);factor=rng.uniform(.1,2)
        def independent(w):
            x,m,v=w;M=b*m+(1-b)*x;V=v+(1-c)*np.tanh(k*(x*x-v))*x*x
            return np.array([x-factor*M/(np.sqrt(V)+e),M,V])
        val,J,det=step_tangent(z,np.eye(3),p,factor)
        cs=np.column_stack([np.imag(independent(z.astype(complex)+1e-25j*np.eye(3)[j]))/1e-25 for j in range(3)])
        errors.append(float(np.max(abs(J-cs))/(1+np.max(abs(cs)))))
        dets.append(float(abs(np.linalg.det(J)-det)/(1+abs(det))))
        values.append(float(np.max(abs(val-independent(z)))/(1+np.max(abs(val)))))
    assert max(errors)<1e-10 and max(dets)<1e-9 and max(values)<1e-13
    return dict(cases=512,complex_step_max_scaled_error=max(errors),determinant_max_scaled_error=max(dets),state_max_scaled_error=max(values))


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    seed=2026110162;u=qmc.Sobol(4,scramble=True,seed=seed).random_base2(12)
    base=np.c_[.02+.975*u[:,0],.02+.975*u[:,1],10**(-6+5*u[:,2]),10**(-3+8*u[:,3])]
    base[2048:,0]=.7+.299*u[2048:,0];base[::2,2]=.1
    extra=[]
    for b,c,e in [(.9806451612903226,.6006451612903226,9.301461332493893e-5),(.9662903225806452,.10709677419354839,9.301461332493893e-5),(.7993178686,.1794898437,1.4848448093e-6),(.9,.999,.1)]:
        for k in [.001,.1,1.,10.,100.,1000.,10000.,100000.]:extra.append([b,c,e,k])
    base=np.r_[base,extra];p=np.repeat(base,2,axis=0);bias=np.tile([False,True],len(base))
    starts=np.tile([.73,0.,0.],(len(p),1));starts[1::2,2]=1e-5/p[1::2,3]
    protocol=dict(seed=seed,declared=len(p),parameter_rows=len(base),objective='theta^2/2',warmup=20000,observed=50000,window=10000,tail=512,
        design='4096 scrambled Sobol rows plus 32 controls. b=.02+.975*u in first half and .7+.299*u in second; c=.02+.975*u; e=.1 on even rows and 10^(-6+5*u) on odd rows; k=10^(-3+8*u).',
        initialization='Each row has x=.73,m=0: even trial v=0, autonomous; odd trial v=1e-5/k, with continuous TFA bias correction sqrt(1-c^t)/(1-b^t). Library initial physical v=1e-6.',
        normalization='alpha=sqrt(k/10), physical epsilon=alpha*e; x=theta/alpha,m=M/alpha,v=V/alpha^2. tanh coefficient in original library is fixed at10.',
        primary_rule='Complete; every leading10k window>.005; tail x span>1e-4; no tail return at lags1..32 under tolerance1e-9*(1+maxabsstate).',
        source='https://github.com/tensorflow/addons/blob/d208d752e98c310280938efa939117bf635a60a8/tensorflow_addons/optimizers/yogi.py',
        scope='Existing library variant, independently transcribed. This scout is not actual TensorFlow execution or a historical-novelty proof.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');(out/'controls.json').write_text(json.dumps(controls(),indent=2)+'\n')
    for name in ['search_yogi_tanh.py','yogi_tanh_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    before=time.perf_counter();records=[]
    for start in range(0,len(p),256):
        stop=min(start+256,len(p));w,t,z,s,x=scout(p[start:stop],starts[start:stop],bias[start:stop],20000,50000,10000,512)
        periods=np.zeros(stop-start,dtype=int)
        for i in range(stop-start):
            if s[i]:continue
            tol=1e-9*(1+np.max(abs(t[i])))
            for lag in range(1,33):
                if np.max(abs(t[i,lag:]-t[i,:-lag]))<tol:periods[i]=lag;break
        block=dict(ids=np.arange(start,stop),parameters=p[start:stop],starts=starts[start:stop],bias=bias[start:stop],windows=w,tails=t,ends=z,status=s,extrema=x,periods=periods)
        np.savez_compressed(out/f'block-{start:05}.npz',**block);records.append(block)
        passed=(s==0)&(w.min(1)>.005)&(np.ptp(t[:,:,0],axis=1)>1e-4)&(periods==0)
        print(stop,'/',len(p),'block flags',int(passed.sum()),'seconds',round(time.perf_counter()-before,1),flush=True)
    d={key:np.concatenate([r[key] for r in records]) for key in records[0]}
    np.savez_compressed(out/'scout.npz',**d)
    good=(d['status']==0)&(d['windows'].min(1)>.005)&(np.ptp(d['tails'][:,:,0],axis=1)>1e-4)&(d['periods']==0)
    result=dict(declared=len(p),completed=int(np.sum(d['status']==0)),primary_ids=np.flatnonzero(good).tolist(),primary_count=int(good.sum()),
        autonomous_flags=int(np.sum(good&~bias)),bias_corrected_flags=int(np.sum(good&bias)),seconds=time.perf_counter()-before)
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())

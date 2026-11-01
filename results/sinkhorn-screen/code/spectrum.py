"""Full tangent spectra and independent Jacobian assembly for selected cases."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import njit
import model

@njit(cache=True)
def state_jacobian(z,L,la,lb,w):
    m,n=L.shape;p=m-1;q=n-1;d=p+q;c=1-w
    x=np.zeros(m);y=np.zeros(n);x[1:]=z[:p];y[1:]=z[p:]
    P=np.empty((m,n));f=np.empty(m)
    for i in range(m):
        s=model.lse(L[i]+y);P[i]=np.exp(L[i]+y-s);f[i]=la[i]-s
    for i in range(1,m):x[i]=c*x[i]+w*(f[i]-f[0])
    Q=np.empty((n,m));g=np.empty(n)
    for j in range(n):
        s=model.lse(L[:,j]+x);Q[j]=np.exp(L[:,j]+x-s);g[j]=lb[j]-s
    for j in range(1,n):y[j]=c*y[j]+w*(g[j]-g[0])
    A=-P[1:,1:]+P[0,1:];B=-Q[1:,1:]+Q[0,1:]
    J=np.zeros((d,d));J[:p,:p]=c*np.eye(p);J[:p,p:]=w*A
    J[p:,:p]=w*c*B;J[p:,p:]=c*np.eye(q)+w*w*(B@A)
    return np.concatenate((x[1:],y[1:])),J

@njit(cache=True)
def spectrum(z,L,w,warm,block,blocks):
    la,lb=model.targets(L);d=len(z);Q=np.eye(d);values=np.zeros((blocks,d))
    for k in range(warm+block*blocks):
        z,J=state_jacobian(z,L,la,lb,w);Q,R=np.linalg.qr(J@Q)
        if k>=warm:values[(k-warm)//block]+=np.log(np.abs(np.diag(R)))/block
    return z,values

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();args.output.mkdir(exist_ok=False,parents=True)
    selected=[(3,43),(4,12),(4,55),(4,64)]
    results=[]
    for n,case in selected:
        r=next(r for r in json.loads((args.input/f'n{n}.json').read_text())['runs'] if r['case']==case)
        L=np.array(r['log_kernel']);z=np.array(r['blocks'][-1]['checkpoint']);w=r['omega'];la,lb=model.targets(L)
        zz,J=state_jacobian(z,L,la,lb,w)
        derivative_error=float(np.max(abs(J-model.jacobian(z,L,la,lb,w))))
        state_error=float(np.max(abs(zz-model.step(z,np.ones(len(z)),L,la,lb,w)[0])))
        assert derivative_error<2e-12 and state_error<2e-12
        t=time.perf_counter();last,rates=spectrum(z,L,w,10000,10000,20)
        expected=len(z)*np.log(abs(1-w));sumerror=float(np.max(abs(rates.sum(1)-expected)))
        assert sumerror<1e-10
        result=dict(case=case,shape=[n,n],omega=w,checkpoint_step=r['blocks'][-1]['stop_step'],starting_state=z.tolist(),
            warm=10000,block_length=10000,blocks=20,rates=rates.tolist(),mean_rates=rates.mean(0).tolist(),
            expected_volume_rate=float(expected),max_block_volume_error=sumerror,
            independent_jacobian_error=derivative_error,state_error=state_error,final_state=last.tolist(),elapsed_seconds=time.perf_counter()-t)
        results.append(result);print(json.dumps(result),flush=True)
    (args.output/'summary.json').write_text(json.dumps(dict(runs=results,sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(model.__file__))}),indent=2)+'\n')

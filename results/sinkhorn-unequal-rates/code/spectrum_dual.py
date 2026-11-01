"""Full QR spectrum with a separately assembled two-block Jacobian."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import njit
import model,model_dual

@njit(cache=True)
def state_jacobian(z,L,la,lb,wx,wy):
    m,n=L.shape;p=m-1;q=n-1;d=p+q;rx=1-wx;ry=1-wy
    x=np.zeros(m);y=np.zeros(n);x[1:]=z[:p];y[1:]=z[p:];P=np.empty((m,n));f=np.empty(m)
    for i in range(m):
        s=model.lse(L[i]+y);P[i]=np.exp(L[i]+y-s);f[i]=la[i]-s
    for i in range(1,m):x[i]=rx*x[i]+wx*(f[i]-f[0])
    Q=np.empty((n,m));g=np.empty(n)
    for j in range(n):
        s=model.lse(L[:,j]+x);Q[j]=np.exp(L[:,j]+x-s);g[j]=lb[j]-s
    for j in range(1,n):y[j]=ry*y[j]+wy*(g[j]-g[0])
    A=-P[1:,1:]+P[0,1:];B=-Q[1:,1:]+Q[0,1:];J=np.zeros((d,d))
    J[:p,:p]=rx*np.eye(p);J[:p,p:]=wx*A;J[p:,:p]=wy*rx*B;J[p:,p:]=ry*np.eye(q)+wx*wy*(B@A)
    return np.concatenate((x[1:],y[1:])),J

@njit(cache=True)
def spectrum(z,L,wx,wy,warm,block,blocks):
    la,lb=model.targets(L);d=len(z);Q=np.eye(d);values=np.zeros((blocks,d))
    for k in range(warm+block*blocks):
        z,J=state_jacobian(z,L,la,lb,wx,wy)
        if not np.isfinite(z).all() or np.max(np.abs(z))>1e5:return z,values,k,False
        Q,R=np.linalg.qr(J@Q)
        if k>=warm:values[(k-warm)//block]+=np.log(np.abs(np.diag(R)))/block
    return z,values,warm+block*blocks,True

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--cases',type=int,nargs='+',required=True)
    args=ap.parse_args();data=json.loads(args.input.read_text());results=[]
    for case in args.cases:
        r=next(r for r in data['runs'] if r['case']==case);L=np.array(r['log_kernel']);z=np.array(r['blocks'][-1]['checkpoint']);wx=r.get('wx',r.get('omega'));wy=r.get('wy',r.get('omega'));la,lb=model.targets(L);d=len(z)
        zz,J=state_jacobian(z,L,la,lb,wx,wy);J2=np.column_stack([model_dual.step(z,v,L,la,lb,wx,wy)[1] for v in np.eye(d)])
        state_error=float(np.max(abs(zz-model_dual.step(z,np.ones(d),L,la,lb,wx,wy)[0])));jacobian_error=float(np.max(abs(J-J2)))
        assert max(state_error,jacobian_error)<2e-12
        t=time.perf_counter();last,rates,steps,complete=spectrum(z,L,wx,wy,10000,10000,20)
        expected=(L.shape[0]-1)*np.log(abs(1-wx))+(L.shape[1]-1)*np.log(abs(1-wy))
        sumerror=float(np.max(abs(rates.sum(1)-expected))) if complete else None
        if complete:assert sumerror<1e-10
        result=dict(case=case,shape=list(L.shape),wx=wx,wy=wy,checkpoint_step=r['blocks'][-1]['stop_step'],starting_state=z.tolist(),warm=10000,block_length=10000,blocks=20,completed=complete,accepted_steps=int(steps),rates=rates.tolist() if complete else None,mean_rates=rates.mean(0).tolist() if complete else None,expected_volume_rate=float(expected),max_block_volume_error=sumerror,independent_jacobian_error=jacobian_error,state_error=state_error,final_state=last.tolist(),elapsed_seconds=time.perf_counter()-t)
        results.append(result);print(json.dumps(result),flush=True)
    args.output.write_text(json.dumps(dict(runs=results,parent_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('model.py'),Path('model_dual.py'))}),indent=2)+'\n')

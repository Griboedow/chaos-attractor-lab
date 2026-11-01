"""Numerical unstable-chart and near-homoclinic proposals, not certificates."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.optimize import brentq
from continue_raar_cycles import iterate
from trace_raar_cycle_seeds import value

def multiply(a,b):return np.convolve(a,b)[:len(a)]

def inverse_sqrt(a):
    n=len(a);b=np.zeros(n);b[0]=np.sqrt(a[0])
    for k in range(1,n):b[k]=(a[k]-sum(b[i]*b[k-i] for i in range(1,k)))/(2*b[0])
    c=np.zeros(n);c[0]=1/b[0]
    for k in range(1,n):c[k]=-sum(b[i]*c[k-i] for i in range(1,k+1))/b[0]
    return c

def projector_jet(x,center,normal,radius):
    v=x.copy();v[:,0]-=center;v-=normal[:,None]*(normal@v)[None,:]
    squared=sum(multiply(row,row) for row in v);reciprocal=inverse_sqrt(squared)
    result=radius*np.array([multiply(row,reciprocal) for row in v]);result[:,0]+=center;return result

def map_jet(x,p):
    beta,r,theta,*c=p;b=projector_jet(x,np.array(c),np.array([0.,-np.sin(theta),np.cos(theta)]),r)
    a=projector_jet(2*b-x,np.zeros(3),np.array([0.,0.,1.]),1.)
    return beta*x+beta*a+(1-2*beta)*b

def chart(center,p,degree=24):
    _,M,_,_,status=iterate(center,p,7);assert status==0
    eigen,basis=np.linalg.eig(M);order=np.argsort(abs(eigen))[::-1];eigen=eigen[order].real;basis=basis[:,order].real
    for i in range(3):
        basis[:,i]/=np.linalg.norm(basis[:,i])
        if basis[np.argmax(abs(basis[:,i])),i]<0:basis[:,i]*=-1
    coeff=np.zeros((3,degree+1));coeff[:,0]=center;coeff[:,1]=basis[:,0]
    for k in range(2,degree+1):
        jet=coeff[:,:k+1].copy()
        for _ in range(7):jet=map_jet(jet,p)
        coeff[:,k]=np.linalg.solve(eigen[0]**k*np.eye(3)-M,jet[:,k])
    return coeff,eigen,basis

@njit(cache=True)
def evaluate(coeff,s):
    x=coeff[:,-1].copy()
    for k in range(coeff.shape[1]-2,-1,-1):x=x*s+coeff[:,k]
    return x

@njit(cache=True,parallel=True)
def trace(coeff,parameters,grid,returns):
    states=np.full((returns+1,len(grid),3),np.nan);status=np.zeros(len(grid),np.int64);axes=np.full(len(grid),np.inf)
    for i in prange(len(grid)):
        x=evaluate(coeff,grid[i]);states[0,i]=x
        for k in range(1,returns+1):
            for j in range(7):
                x,flag,axis=value(x,parameters);axes[i]=min(axes[i],axis)
                if flag:status[i]=flag;break
            if status[i]:break
            states[k,i]=x
    return states,status,axes

def main(args):
    repo=Path(__file__).resolve().parents[1];out=args.out;out.mkdir(parents=True,exist_ok=False);begin=time.time();set_num_threads(8)
    certpath=repo/'results/raar-cycle-structure/cycle-proof/base-unstable-period-7.json';cert=json.loads(certpath.read_text());p=np.array(cert['parameters']);center=np.array([float(v) for v in cert['center_decimal']]);coeff,eigen,basis=chart(center,p);inverse=np.linalg.inv(basis)
    residual=[]
    for s in np.linspace(-.02,.02,101):residual.append(float(np.linalg.norm(iterate(evaluate(coeff,s),p,7)[0]-evaluate(coeff,eigen[0]*s))))
    assert max(residual)<1e-11
    grid=np.linspace(-.02,.02,8193);states,status,axes=trace(coeff,p,grid,32);assert not status.any();coordinates=(states-center)@inverse.T
    np.savez_compressed(out/'trace.npz',params=p,center=center,coefficients=coeff,multipliers=eigen,basis=basis,grid=grid,states=states,coordinates=coordinates,status=status,minimum_axes=axes)
    attempts=[]
    for k in range(3,33):
        uv=coordinates[k,:,0];pairs=np.where(uv[:-1]*uv[1:]<0)[0]
        for i in pairs:
            if min(np.linalg.norm(coordinates[k,i,1:]),np.linalg.norm(coordinates[k,i+1,1:]))>.1:continue
            if max(np.linalg.norm(states[:k+1,i]-center,axis=1))<.1:continue
            if np.linalg.norm(states[k,i+1]-states[k,i])>.2:continue
            try:
                fun=lambda s:float((inverse@(iterate(evaluate(coeff,s),p,7*k)[0]-center))[0])
                t=brentq(fun,grid[i],grid[i+1],xtol=1e-16,rtol=1e-14);end,M,gap,orbit,flag=iterate(evaluate(coeff,t),p,7*k);coord=inverse@(end-center)
                row=dict(returns=k,grid_index=int(i),parameter=float(t),end_coordinates=coord.tolist(),end_distance=float(np.linalg.norm(end-center)),initial_distance=float(np.linalg.norm(evaluate(coeff,t)-center)),max_excursion=float(np.max(np.linalg.norm(orbit-center,axis=1))),minimum_axis=float(gap),root_residual=float(abs(coord[0])),status=int(flag))
                attempts.append(row)
            except Exception as e:attempts.append(dict(returns=k,grid_index=int(i),error=repr(e)))
        good=[r for r in attempts if r.get('returns')==k and 'end_distance' in r]
        if good:print(json.dumps(dict(returns=k,roots=len(good),closest=min(good,key=lambda r:r['end_distance']))),flush=True)
    ordered=sorted([r for r in attempts if 'end_distance' in r],key=lambda r:r['end_distance'])
    report=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cycle_certificate_sha256=hashlib.sha256(certpath.read_bytes()).hexdigest(),max_chart_residual=max(residual),chart_degree=24,multipliers=eigen.tolist(),basis_condition=float(np.linalg.cond(basis)),attempts=attempts,closest=ordered[:30],seconds=time.time()-begin,scope='Numerical proposal only. The polynomial is not a validated invariant chart; the final linear stable plane is not the stable manifold. No homoclinic connection or chaos theorem.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps(dict(total_attempts=len(attempts),best=ordered[:3],seconds=report['seconds'])),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);main(parser.parse_args())

"""Continue all seven rotating-frame sectors of the nonsymmetric 2T orbit by pseudo-arclength."""
import argparse
import json
import os
import sys
from pathlib import Path

for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'flows'))
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from adaptive import polar_ring
from field import ring_tangent
from relative_cores import rotation

E0=-.4297716347835525
DELTA=.0759375
R=rotation(2*np.pi*4/7)


@njit(cache=True)
def variational(t,y,p):
    out=np.empty_like(y)
    for j in range(len(y)//15):
        k=j*15
        v=np.empty(6);v[:3]=y[k:k+3]
        for c in range(4):
            v[3:]=y[k+3+c*3:k+6+c*3];d=ring_tangent(v,p)
            out[k+3+c*3:k+6+c*3]=d[3:]
        out[k:k+3]=d[:3];ss=v[0]*v[0]+v[1]*v[1]
        out[k+14]+=DELTA*(ss-1)/(ss+1)
    return out


def equations(x,p0,phase,level,tol=2e-10):
    m=(len(x)-2)//3;u=x[:3*m].reshape(m,3);tau,f=x[-2:];p=p0.copy();p[8]=E0+f*DELTA
    y=np.zeros((m,15));y[:,:3]=u;y[:,3:12]=np.eye(3).reshape(9)
    sol=solve_ivp(variational,(0,tau),y.ravel(),args=(p,),method='DOP853',rtol=tol,atol=tol*.01)
    if not sol.success:raise RuntimeError(sol.message)
    end=sol.y[:,-1].reshape(m,15);res=np.zeros(3*m+1);J=np.zeros((3*m+1,len(x)));M=np.eye(3)
    for i,row in enumerate(end):
        indices=slice(3*i,3*i+3);j=(i+1)%m;nxt=slice(3*j,3*j+3)
        A=R.T@row[3:12].reshape(3,3).T
        res[indices]=R.T@row[:3]-u[j];J[indices,indices]+=A;J[indices,nxt]-=np.eye(3)
        J[indices,-2]=R.T@polar_ring(row[:3],p);J[indices,-1]=R.T@row[12:15];M=A@M
    res[-1]=x[phase]-level;J[-1,phase]=1
    if m==1:M=np.linalg.matrix_power(M,7)
    eig=np.linalg.eigvals(M);neutral=np.argmin(abs(eig-1));other=np.delete(eig,neutral)
    mu=other[np.argmax(abs(other))]
    return res,J,[float(mu.real),float(mu.imag)],float(eig[neutral].real)


def solve(x,p,phase,level,normal,point,tol):
    for _ in range(10):
        res,J,mu,one=equations(x,p,phase,level,tol)
        full=np.r_[res,np.dot(normal,x-point)]
        if max(abs(full))<tol*.8:return x,J,mu,one,float(max(abs(res)))
        mat=np.vstack([J,normal])
        try:step=np.linalg.solve(mat,-full)
        except np.linalg.LinAlgError:return None
        if np.linalg.norm(step)>.10:return None
        x=x+step
    return None


def continuation(case,p,z,T,phase,direction,steps,ds,tol,path):
    m=7 if case=='asymmetric' else 1
    locs=[]
    init=solve_ivp(lambda t,s:polar_ring(s,p),(0,T),z,method='DOP853',rtol=3e-13,atol=3e-15,dense_output=True)
    for j in range(m):locs.append(rotation(-2*np.pi*4*j/7)@init.sol(j*T/7))
    x=np.r_[np.concatenate(locs),T/7,(p[8]-E0)/DELTA];normal=np.zeros(len(x));normal[-1]=1
    got=solve(x,p,phase,z[phase],normal,x,tol)
    if got is None:raise RuntimeError('Initial correction failed')
    x,J,mu,one,res=got;v=np.linalg.svd(J,full_matrices=True)[2][-1]
    if v[-1]*direction<0:v=-v
    initial=x.copy();total=0.;h=ds
    with path.open('w') as file:
        for it in range(steps+1):
            u=x[:-2].reshape(m,3);spread=np.linalg.norm(u-u.mean(0),axis=1).max()
            row=dict(i=it,branch=case,arclength=total,parameter_fraction=x[-1],e=E0+DELTA*x[-1],
                     period=7*x[-2],multiplier=mu,neutral=one,relative_node_spread=spread,
                     de_ds=DELTA*v[-1],residual=res,phase=phase,phase_level=z[phase],nodes=u.tolist())
            file.write(json.dumps(row,separators=(',',':'))+'\n');file.flush()
            if it%25==0:print(case,it,x[-1],mu[0],spread,flush=True)
            if it==steps or (it>20 and (x[-1]<.86 or x[-1]>1.17)):break
            while h>=1e-6:
                pred=x+h*v
                got=solve(pred,p,phase,z[phase],v,pred,tol)
                if got is not None:break
                h*=.5
            if got is None:break
            prev=x;x,J,mu,one,res=got;w=np.linalg.svd(J,full_matrices=True)[2][-1]
            if np.dot(w,v)<0:w=-w
            total+=np.linalg.norm(x-prev);curve=np.dot(w,v);v=w
            h=min(ds*1.4,h*(1.20 if curve>.97 else .65))
            if it>20 and total>.14 and np.linalg.norm(x-initial)<3e-4:break


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--old',required=True,type=Path);ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--branch',choices=['asymmetric','symmetric'],required=True)
    ap.add_argument('--direction',type=int,choices=[-1,1],default=-1)
    ap.add_argument('--steps',type=int,default=310);ap.add_argument('--ds',type=float,default=.00075)
    ap.add_argument('--tol',type=float,default=4e-10)
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    if a.branch=='asymmetric':q=next(r for r in json.loads((a.old/'nested-skeletons.json').read_text())['cores'] if r['period_order']==2)
    else:q=json.loads((a.old/'final-link-crosscheck.json').read_text())['symmetric_companion']
    z=np.array(q['initial']);p=np.array(q['parameters']);phase=int(np.argmax(abs(polar_ring(z,p))))
    continuation(a.branch,p,z,q['period'],phase,a.direction,a.steps,a.ds,a.tol,a.out/'branch.jsonl')


if __name__=='__main__':
    main()

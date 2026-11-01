"""Full variational equations for the proper-tetrahedral vortex flow."""
from pathlib import Path
import sys,json,time
import numpy as np
from numba import njit,prange,set_num_threads
sys.path.insert(0,str(Path('D:/Codex/chaos-attractor-lab/research')))
from twist_map import orthonormalize
from tetra_vortex import rhs
set_num_threads(4)

@njit
def field(z,centers,axes,sigma,gamma):
    f=-gamma*z;jac=-gamma*np.eye(3);sig2=sigma*sigma
    for i in range(len(centers)):
        r=z-centers[i];u=axes[i]
        a=u[1]*r[2]-u[2]*r[1];b=u[2]*r[0]-u[0]*r[2];c=u[0]*r[1]-u[1]*r[0]
        g=np.exp(-.5*np.dot(r,r)/sig2);cross=np.array([a,b,c]);f+=g*cross
        for j in range(3):
            for k in range(3):jac[j,k]-=g*cross[j]*r[k]/sig2
        jac[0,1]-=g*u[2];jac[0,2]+=g*u[1];jac[1,0]+=g*u[2];jac[1,2]-=g*u[0];jac[2,0]-=g*u[1];jac[2,1]+=g*u[0]
    return f,jac

@njit
def advance(z,q,c,u,s,g,h):
    a,j=field(z,c,u,s,g);av=j@q
    b,j=field(z+.5*h*a,c,u,s,g);bv=j@(q+.5*h*av)
    d,j=field(z+.5*h*b,c,u,s,g);dv=j@(q+.5*h*bv)
    e,j=field(z+h*d,c,u,s,g);ev=j@(q+h*dv)
    return z+h/6*(a+2*b+2*d+e),q+h/6*(av+2*bv+2*dv+ev)

@njit(parallel=True)
def integrate(centers,axes,sigma,gamma,starts,hs,warm_time,observed_time,window_time,keep):
    count=len(starts);w=np.zeros((count,int(observed_time/window_time),3));tails=np.zeros((count,keep,3));means=np.zeros_like(w);maxradius=np.zeros(count);ends=np.zeros_like(starts)
    for i in prange(count):
        h=hs[i];warm=int(round(warm_time/h));observed=int(round(observed_time/h));window=int(round(window_time/h));stride=observed//keep
        z=starts[i].copy();q=np.eye(3);acc=np.zeros(3);mean=np.zeros(3)
        for t in range(warm+observed):
            z,q=advance(z,q,centers[i],axes[i],sigma[i],gamma[i],h);q,logs=orthonormalize(q)
            if t>=warm:
                acc+=logs;mean+=z;maxradius[i]=max(maxradius[i],np.linalg.norm(z))
                if (t-warm+1)%window==0:
                    ii=(t-warm+1)//window-1;w[i,ii]=acc/(window*h);means[i,ii]=mean/window;acc[:]=0;mean[:]=0
                if (t-warm+1)%stride==0:tails[i,(t-warm+1)//stride-1]=z
        ends[i]=z
    return w,tails,means,maxradius,ends

def main():
    root=Path(__file__).parent;d=np.load(root/'tetra-v1/data.npz');out=root/'tetra-spectra';out.mkdir(exist_ok=False)
    ids=np.repeat([421,565,705,117],3);hs=np.tile([.05,.025,.0125],4)
    ferr=0.;jerr=0.;diverr=0.;rng=np.random.default_rng(2026110136)
    for i in [421,565,705,117]:
        for z in rng.normal(size=(20,3)):
            f,j=field(z,d['centers'][i],d['axes'][i],d['sigma'][i],d['gamma'][i]);ff,jv=rhs(z,np.ones(3),d['centers'][i],d['axes'][i],d['sigma'][i],d['gamma'][i])
            ferr=max(ferr,float(np.max(np.abs(f-ff))));jerr=max(jerr,float(np.max(np.abs(j@np.ones(3)-jv))));diverr=max(diverr,abs(float(np.trace(j)+3*d['gamma'][i])))
            eps=1e-5
            jj=np.column_stack([(field(z+eps*e,d['centers'][i],d['axes'][i],d['sigma'][i],d['gamma'][i])[0]-field(z-eps*e,d['centers'][i],d['axes'][i],d['sigma'][i],d['gamma'][i])[0])/(2*eps) for e in np.eye(3)])
            jerr=max(jerr,float(np.max(np.abs(j-jj))))
    assert ferr<1e-13 and jerr<1e-8 and diverr<1e-13
    checks={'field_error':ferr,'jacobian_error':jerr,'divergence_error':diverr};(out/'checks.json').write_text(json.dumps(checks,indent=2))
    begin=time.time();w,z,means,bounds,ends=integrate(d['centers'][ids],d['axes'][ids],d['sigma'][ids],d['gamma'][ids],d['starts'][ids],hs,3000,10000,1000,100000)
    np.savez_compressed(out/'data.npz',ids=ids,h=hs,windows=w,tail=z,means=means,max_radius=bounds,end=ends)
    rows=[{'id':int(i),'h':float(h),'spectrum':a.mean(0).tolist(),'min_leading':float(a[:,0].min()),'sum_error':float(np.max(np.abs(a.sum(1)+3*d['gamma'][i]))),'radius':float(r)} for i,h,a,r in zip(ids,hs,w,bounds)]
    (out/'summary.json').write_text(json.dumps({'seconds':time.time()-begin,'rows':rows},indent=2));print(json.dumps(rows),flush=True)

if __name__=='__main__':main()

from paths import OUTPUT, ROOT, RESEARCH
"""Tetrahedrally equivariant localized vortices with linear confinement."""
from pathlib import Path
import itertools,json,time
import numpy as np
from scipy.stats import qmc
from numba import njit,prange,set_num_threads
set_num_threads(4)

def tetrahedral_group():
    group=[]
    for shift in range(3):
        cyclic=np.roll(np.eye(3),shift,axis=0)
        for signs in itertools.product([-1.,1.],repeat=3):
            if np.prod(signs)>0:group.append(np.diag(np.array(signs))@cyclic)
    return np.array(group)

@njit
def rhs(z,v,centers,axes,sigma,gamma):
    f=-gamma*z;dv=-gamma*v;sig2=sigma*sigma
    for j in range(len(centers)):
        r=z-centers[j];u=axes[j]
        cr=np.cross(u,r);cv=np.cross(u,v)
        g=np.exp(-.5*np.dot(r,r)/sig2)
        f+=g*cr;dv+=g*(cv-cr*np.dot(r,v)/sig2)
    return f,dv

@njit
def advance(z,v,c,u,s,g,h):
    a,av=rhs(z,v,c,u,s,g);b,bv=rhs(z+h*.5*a,v+h*.5*av,c,u,s,g)
    d,dv=rhs(z+h*.5*b,v+h*.5*bv,c,u,s,g);e,ev=rhs(z+h*d,v+h*dv,c,u,s,g)
    return z+h/6*(a+2*b+2*d+e),v+h/6*(av+2*bv+2*dv+ev)

@njit(parallel=True)
def survey(centers,axes,sigma,gamma,starts,warm,observed,h,keep):
    count=len(starts);windows=np.zeros((count,4));tails=np.zeros((count,keep,3));ends=np.zeros((count,3));speed=np.zeros(count);radius=np.zeros(count)
    for lane in prange(count):
        z=starts[lane].copy();v=np.array([1.,0.,0.]);acc=0.
        for t in range(warm+observed):
            z,v=advance(z,v,centers[lane],axes[lane],sigma[lane],gamma[lane],h)
            norm=np.linalg.norm(v);v/=norm
            if t>=warm:
                radius[lane]=max(radius[lane],np.linalg.norm(z));acc+=np.log(norm)
                if (t-warm+1)%(observed//4)==0:
                    windows[lane,(t-warm+1)//(observed//4)-1]=acc/(h*(observed//4));acc=0.
                if (t-warm+1)%(observed//keep)==0:tails[lane,(t-warm+1)//(observed//keep)-1]=z
        ends[lane]=z;speed[lane]=np.linalg.norm(rhs(z,v,centers[lane],axes[lane],sigma[lane],gamma[lane])[0])
    return windows,tails,ends,speed,radius

def main():
    out=OUTPUT/'tetra-v1';out.mkdir(exist_ok=False)
    group=tetrahedral_group();u=qmc.Sobol(9,scramble=True,seed=2026110132).random_base2(10)
    seed_c=np.column_stack([np.ones(len(u)),.05+.85*u[:,0],.05+.85*u[:,1]])
    seed_c/=np.linalg.norm(seed_c,axis=1)[:,None]
    phi=2*np.pi*u[:,4];cost=2*u[:,5]-1
    seed_axis=np.column_stack([np.sqrt(1-cost**2)*np.cos(phi),np.sqrt(1-cost**2)*np.sin(phi),cost])
    centers=np.einsum('gij,nj->ngi',group,seed_c);axes=np.einsum('gij,nj->ngi',group,seed_axis)
    sigma=.18+.75*u[:,2];gamma=10**(-3+2.3*u[:,3]);starts=(2*u[:,6:9]-1)*1.2
    # Algebra checks precede the numerical survey.
    residual=0.
    for i in range(8):
        z=starts[i];v=np.array([1.,0.,0.]);f=rhs(z,v,centers[i],axes[i],sigma[i],gamma[i])[0]
        for r in group:residual=max(residual,float(np.max(np.abs(rhs(r@z,v,centers[i],axes[i],sigma[i],gamma[i])[0]-r@f))))
    assert residual<1e-13
    begin=time.time();w,tail,end,speed,radius=survey(centers,axes,sigma,gamma,starts,4000,16000,.05,4000)
    np.savez_compressed(out/'data.npz',group=group,centers=centers,axes=axes,sigma=sigma,gamma=gamma,starts=starts,windows=w,tail=tail,end=end,speed=speed,radius=radius)
    ids=np.where((w.min(1)>.003)&(speed>1e-4))[0]
    print(json.dumps({'seconds':time.time()-begin,'equivariance_error':residual,'positive':ids.tolist(),'max_leading':float(w.mean(1).max()),'nontrivial_speed':int((speed>1e-4).sum())}),flush=True)

if __name__=='__main__':main()

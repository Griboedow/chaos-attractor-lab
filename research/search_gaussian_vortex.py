from pathlib import Path
import argparse,json,time
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.stats import qmc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
set_num_threads(4)

# Three orthogonal localized vortex axes; centers are off their axes through zero.
CENTERS=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
AXES=np.array([[0.,1.,0.],[0.,0.,1.],[1.,0.,0.]])

@njit
def rhs(z,v,p):
    f=-p[0]*z+p[8:11]
    dv=-p[0]*v
    for n in range(3):
        rr=z-p[7]*CENTERS[n]
        u=AXES[n]
        cross=np.cross(u,rr);cv=np.cross(u,v)
        g=p[1+n]*np.exp(-.5*np.dot(rr,rr)/p[4+n]**2)
        f+=g*cross
        dv+=g*(cv-cross*np.dot(rr,v)/p[4+n]**2)
    return f,dv

@njit
def advance(z,v,p,h):
    a,av=rhs(z,v,p)
    b,bv=rhs(z+h*.5*a,v+h*.5*av,p)
    c,cv=rhs(z+h*.5*b,v+h*.5*bv,p)
    d,dv=rhs(z+h*c,v+h*cv,p)
    return z+h/6*(a+2*b+2*c+d),v+h/6*(av+2*bv+2*cv+dv)

@njit(parallel=True)
def survey(parameters,warm,observed,h,keep):
    count=len(parameters);w=np.zeros((count,4));tail=np.zeros((count,keep,3));end=np.zeros((count,3));speed=np.zeros(count)
    for i in prange(count):
        z=np.zeros(3);z[0]=.1;z[1]=.2;z[2]=.3
        v=np.zeros(3);v[0]=1
        accum=0.
        stride=observed//keep
        for t in range(warm+observed):
            z,v=advance(z,v,parameters[i],h)
            norm=np.sqrt(np.dot(v,v));v/=norm
            if t>=warm:
                accum+=np.log(norm)
                if (t-warm+1)%(observed//4)==0:
                    w[i,(t-warm+1)//(observed//4)-1]=accum/(observed//4*h);accum=0.
                if (t-warm+1)%stride==0:
                    tail[i,(t-warm+1)//stride-1]=z
        end[i]=z
        f,d=rhs(z,v,parameters[i]);speed[i]=np.linalg.norm(f)
    return w,tail,end,speed

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True)
    out=parser.parse_args().out;out.mkdir(parents=True,exist_ok=False)
    u=qmc.Sobol(11,scramble=True,seed=20261108).random_base2(10)
    p=np.empty_like(u)
    p[:,0]=.03+.6*u[:,0]
    p[:,1:4]=(u[:,1:4]*2-1)*12
    p[:,4:7]=.45+1.0*u[:,4:7]
    p[:,7]=.6+1.4*u[:,7]
    p[:,8:11]=(u[:,8:11]*2-1)*.7
    # A fully symmetric quarter, followed by anisotropic widths/strengths.
    p[:256,2:4]=p[:256,1,None];p[:256,5:7]=p[:256,4,None];p[:256,9:11]=p[:256,8,None]
    t=time.time();w,z,end,speed=survey(p,10000,20000,.025,2000)
    np.savez_compressed(out/'data.npz',parameters=p,windows=w,tail=z,end=end,speed=speed)
    ids=np.where((w.min(1)>.005)&(speed>1e-4))[0]
    print(json.dumps({'seconds':time.time()-t,'positive':len(ids),'ids':ids.tolist()}),flush=True)
    for page in range((len(ids)+15)//16):
        fig,axes=plt.subplots(4,4,figsize=(12,12),facecolor='#040813')
        for ax,i in zip(axes.flat,ids[page*16:page*16+16]):
            vals,basis=np.linalg.eigh(np.cov(z[i].T));points=(z[i]-z[i].mean(0))@basis[:,::-1]
            ax.scatter(points[:,0],points[:,1],s=.18,c=points[:,2],cmap='magma')
            ax.set_aspect('equal');ax.set_axis_off();ax.set_facecolor('black');ax.set_title(f'{i} L={w[i].mean():.3f}',color='white')
        fig.tight_layout();fig.savefig(out/f'gallery-{page}.png',dpi=120);plt.close(fig)

if __name__=='__main__':main()

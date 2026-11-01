import os
os.environ.setdefault('MPLCONFIGDIR','D:/Codex/twist-search-20261101/mpl')
import argparse,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.stats import qmc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
set_num_threads(4)

@njit
def advance(z,v,p):
 r=p[0]
 z=r*z;z[0]+=1-r;v=r*v
 for stage in range(3):
  i=stage;j=(stage+1)%3;k=(stage+2)%3
  w=z[k];angle=p[1+stage]+p[4+stage]*w+p[7+stage]*w*w*w
  derivative=p[4+stage]+3*p[7+stage]*w*w
  c=np.cos(angle);s=np.sin(angle)
  u=c*z[i]-s*z[j];vv=s*z[i]+c*z[j]
  du=c*v[i]-s*v[j]-vv*derivative*v[k]
  dv=s*v[i]+c*v[j]+u*derivative*v[k]
  z[i]=u;z[j]=vv;v[i]=du;v[j]=dv
 return z,v

@njit(parallel=True)
def survey(parameters,warmup,observed,keep):
 n=len(parameters);tail=np.zeros((n,keep,3));windows=np.zeros((n,4));ends=np.zeros((n,3))
 for i in prange(n):
  z=np.zeros(3);z[0]=.1;z[1]=.2;z[2]=.3
  v=np.zeros(3);v[0]=1.;v[1]=.1;v[2]=.2
  v/=np.sqrt((v*v).sum());acc=0.
  for t in range(warmup+observed):
   z,v=advance(z,v,parameters[i]);norm=np.sqrt((v*v).sum());v/=norm
   if t>=warmup:
    acc+=np.log(norm)
    if (t-warmup+1)%(observed//4)==0:
     windows[i,(t-warmup+1)//(observed//4)-1]=acc/(observed//4);acc=0.
    if t>=warmup+observed-keep:tail[i,t-(warmup+observed-keep)]=z
  ends[i]=z
 return tail,windows,ends

def features(z):
 vals=np.linalg.eigvalsh(np.cov(z.T));center=z.mean(0);scale=np.ptp(z,axis=0)
 if np.min(scale)<1e-10:return np.zeros(13)
 zz=(z-z.min(0))/scale
 occ=[]
 for grid in [8,16,32]:
  bins=np.minimum((zz*grid).astype(int),grid-1);occ.append(len(np.unique(bins,axis=0)))
 return np.r_[vals/(vals.sum()+1e-30),center,scale,np.array(occ)/len(z),np.mean(np.linalg.norm(z-center,axis=1))]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--power',type=int,default=12);args=ap.parse_args()
 args.output.mkdir(parents=True,exist_ok=False);t=time.time()
 u=qmc.Sobol(10,scramble=True,seed=20261104).random_base2(args.power)
 p=np.empty_like(u);p[:,0]=.35+.58*u[:,0];p[:,1:4]=2*np.pi*u[:,1:4]-np.pi;p[:,4:7]=24*u[:,4:7]-12;p[:,7:]=16*u[:,7:]-8
 # First half uses only linear-coordinate twists, second half permits cubics.
 p[:len(p)//2,7:]=0
 tail,w,end=survey(p,5000,12000,2048)
 f=np.array([features(z) for z in tail]);valid=(w.min(1)>.008)&(w.mean(1)<.55)&(f[:,0]>.004)&(f[:,9]<.45)
 ids=np.where(valid)[0]
 np.savez_compressed(args.output/'pilot.npz',parameters=p,windows=w,end=end,features=f,valid_ids=ids,tail=tail)
 # Farthest-point selection of raw moment/occupancy summaries encourages variety.
 ff=f[ids];ff=(ff-ff.mean(0))/(ff.std(0)+1e-8)
 selected=[int(np.argmin(w[ids].mean(1)))];distance=np.full(len(ids),np.inf)
 for k in range(min(64,len(ids))-1):
  distance=np.minimum(distance,np.linalg.norm(ff-ff[selected[-1]],axis=1));distance[selected]=-1
  selected.append(int(np.argmax(distance)))
 chosen=ids[selected];(args.output/'selection.json').write_text(json.dumps(chosen.tolist()))
 for page in range((len(chosen)+15)//16):
  fig=plt.figure(figsize=(14,14),facecolor='#080f1e')
  for j,i in enumerate(chosen[page*16:page*16+16]):
   ax=fig.add_subplot(4,4,j+1,projection='3d',facecolor='#080f1e');z=tail[i]
   ax.scatter(*z.T,c=np.arange(len(z)),cmap='turbo',s=.4,alpha=.7)
   ax.set_title(f'{i}  L={w[i].mean():.3f}  r={p[i,0]:.2f}',color='white',fontsize=10)
   ax.set_axis_off();ax.view_init(23,-58);ax.set_box_aspect(np.ptp(z,axis=0))
  fig.tight_layout();fig.savefig(args.output/f'gallery-{page}.png',dpi=130);plt.close(fig)
 print(json.dumps({'lanes':len(p),'positive':int((w.min(1)>.008).sum()),'geometry_candidates':len(ids),'selected':chosen.tolist(),'seconds':time.time()-t}),flush=True)

if __name__=='__main__':main()

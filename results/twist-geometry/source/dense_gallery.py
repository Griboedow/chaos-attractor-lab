import os
os.environ.setdefault('MPLCONFIGDIR','D:/Codex/twist-search-20261101/mpl')
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pilot import survey
base=Path(__file__).parent
source=np.load(base/'pilot-v2/pilot.npz')
p=source['parameters'];w=source['windows'];f=source['features']
valid=(w.min(1)>.025)&(w.mean(1)<1.15)&(f[:,0]>.007)&(f[:,10]<.55)
ids=np.where(valid)[0];ff=f[ids];ff=(ff-ff.mean(0))/(ff.std(0)+1e-8)
chosen=[int(np.argmin(w[ids].mean(1)))];dist=np.full(len(ids),np.inf)
for k in range(95):
 dist=np.minimum(dist,np.linalg.norm(ff-ff[chosen[-1]],axis=1));dist[chosen]=-1;chosen.append(int(np.argmax(dist)))
ids=ids[chosen]
out=base/'dense';out.mkdir(exist_ok=False)
z,windows,end=survey(p[ids],10000,60000,50000)
np.savez_compressed(out/'data.npz',ids=ids,parameters=p[ids],tail=z,windows=windows,end=end)
for page in range(6):
 fig,axes=plt.subplots(4,4,figsize=(13,13),facecolor='#050914')
 for j,ax in enumerate(axes.flat):
  idx=page*16+j;points=z[idx];center=points.mean(0)
  vals,vec=np.linalg.eigh(np.cov(points.T));coords=(points-center)@vec[:,::-1]
  # Orthographic view after a rigid PCA rotation; never a nonlinear embedding.
  az=.4;el=.35;rot=np.array([[np.cos(az),0,np.sin(az)],[np.sin(az)*np.sin(el),np.cos(el),-np.cos(az)*np.sin(el)],[0,0,1]])
  xy=(coords@rot.T)[:,:2];lim=np.max(np.ptp(xy,axis=0))*.52;mid=(xy.max(0)+xy.min(0))/2
  h,xe,ye=np.histogram2d(xy[:,0]-mid[0],xy[:,1]-mid[1],bins=350,range=[[-lim,lim],[-lim,lim]])
  ax.imshow(np.log1p(h.T),origin='lower',cmap='magma',interpolation='nearest',vmin=0)
  ax.set_title(f'{ids[idx]}  L={windows[idx].mean():.3f}  r={p[ids[idx],0]:.2f}',color='white',fontsize=10);ax.set_axis_off()
 fig.tight_layout();fig.savefig(out/f'gallery-{page}.png',dpi=150);plt.close(fig)
print('selected',ids.tolist(),flush=True)

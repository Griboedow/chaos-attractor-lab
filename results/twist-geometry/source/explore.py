from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pilot import survey

base=Path(__file__).parent
src=np.load(base/'pilot-v2/pilot.npz')
ids=np.array([307,3254,1709,3555,815,2259])
p=src['parameters'][ids]
out=base/'views';out.mkdir(exist_ok=True)
z,w,end=survey(p,100000,1000000,800000)
np.savez_compressed(out/'data.npz',ids=ids,parameters=p,tail=z,windows=w,end=end)
for ii,ident in enumerate(ids):
 points=z[ii];center=points.mean(0);vals,vec=np.linalg.eigh(np.cov(points.T));basis=vec[:,::-1]
 if np.linalg.det(basis)<0:basis[:,-1]*=-1
 coords=(points-center)@basis
 fig,axes=plt.subplots(2,3,figsize=(15,10),facecolor='#03050b')
 for ax,az in zip(axes.flat,np.arange(6)*np.pi/3):
  el=.35
  projection=np.array([[np.cos(az),0,np.sin(az)],[np.sin(az)*np.sin(el),np.cos(el),-np.cos(az)*np.sin(el)]])
  xy=coords@projection.T;mid=(xy.max(0)+xy.min(0))/2;lim=np.max(np.ptp(xy,axis=0))*.54
  h,*_=np.histogram2d(*(xy-mid).T,bins=700,range=[[-lim,lim]]*2)
  ax.imshow(np.log1p(h.T),origin='lower',cmap='magma',interpolation='nearest',vmin=0)
  ax.set_axis_off();ax.set_title(f'{ident} / {az*180/np.pi:.0f} deg / L={w[ii].mean():.4f}',color='white')
 fig.tight_layout();fig.savefig(out/f'{ident}.png',dpi=125);plt.close(fig)
print(json.dumps({'ids':ids.tolist(),'windows':w.tolist()}),flush=True)

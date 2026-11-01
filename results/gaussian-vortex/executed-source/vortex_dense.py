from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from vortex_pilot import survey
base=Path(__file__).parent
src=np.load(base/'vortex-v1/data.npz')
ids=np.where((src['windows'].min(1)>.005)&(src['speed']>1e-4))[0]
out=base/'vortex-dense';out.mkdir(exist_ok=False)
w,z,end,speed=survey(src['parameters'][ids],40000,160000,.0125,40000)
np.savez_compressed(out/'data.npz',ids=ids,parameters=src['parameters'][ids],windows=w,tail=z,end=end,speed=speed)
for page in range((len(ids)+15)//16):
 fig,axes=plt.subplots(4,4,figsize=(13,13),facecolor='#040813')
 for ax in axes.flat:ax.set_axis_off();ax.set_facecolor('black')
 for ax,j in zip(axes.flat,range(page*16,min((page+1)*16,len(ids)))):
  v,b=np.linalg.eigh(np.cov(z[j].T));points=(z[j]-z[j].mean(0))@b[:,::-1]
  a=.4;e=.35;proj=np.array([[np.cos(a),0,np.sin(a)],[np.sin(a)*np.sin(e),np.cos(e),-np.cos(a)*np.sin(e)]])
  xy=points@proj.T;xy-=(xy.min(0)+xy.max(0))/2;lim=np.max(np.ptp(xy,axis=0))*.53
  h,*_=np.histogram2d(*xy.T,bins=400,range=[[-lim,lim]]*2)
  ax.imshow(np.log1p(h.T),origin='lower',cmap='magma');ax.set_title(f'{ids[j]} L={w[j].mean():.3f}',color='white')
 fig.tight_layout();fig.savefig(out/f'gallery-{page}.png',dpi=125);plt.close(fig)
print('positive',ids[w.min(1)>.005].tolist(),flush=True)

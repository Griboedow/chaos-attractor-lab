import sys,json,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')
from adabelief_map import spectrum
pairs=np.array([[.9,.99],[.9,.999],[.8,.99],[.5,.9],[.5,.13],[.52,.15],
 [.4,.58],[.43,.64],[.32,.05],[.9,.9],[.95,.98],[.8,.8]])
z=np.tile([.73,0.,0.],(len(pairs),1))
t=time.time()
r=spectrum(z,*pairs.T,steps=100000,keep=25000)
out=Path(__file__).parent
np.savez_compressed(out/'gallery.npz',parameters=pairs,**r)
print(json.dumps({'seconds':time.time()-t,'spectra':r['spectrum'].tolist(),'volume_error':float(abs(r['spectrum'].sum(axis=1)-r['determinant_sum']).max())},indent=2),flush=True)
fig=plt.figure(figsize=(16,12),facecolor='#0b1120')
for i,p in enumerate(pairs):
 ax=fig.add_subplot(3,4,i+1,projection='3d',facecolor='#0b1120')
 x,m,s=r['tail'][:,i,:].T
 ax.scatter(x,m,s,c=np.arange(len(x)),s=.12,cmap='turbo',alpha=.5)
 ax.set_title(f'b={p[0]} c={p[1]} | L={r["spectrum"][i,0]:.4f}',color='white',fontsize=10)
 ax.set_xlabel('x',color='white');ax.set_ylabel('m',color='white');ax.set_zlabel('s',color='white')
 ax.tick_params(colors='#8b95a7',labelsize=6)
 for a in [ax.xaxis,ax.yaxis,ax.zaxis]: a.pane.fill=False
fig.tight_layout();fig.savefig(out/'gallery.png',dpi=140)

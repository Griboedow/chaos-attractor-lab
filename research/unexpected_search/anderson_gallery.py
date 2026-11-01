from paths import OUTPUT, ROOT, RESEARCH
"""Actual three-lag states; no nonlinear geometric transform."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=OUTPUT
d=np.load(root/'anderson-v2/data.npz')
ids=d['positive_ids'];tails=d['positive_tails'];rows=[]
for i,z in zip(ids,tails):
    a=z-z.mean(0);cov=np.linalg.eigvalsh(np.cov(a.T));rows.append({'id':int(i),'parameters':d['parameters'][i].tolist(),'L':float(d['windows'][i].mean()),'cov_min_max':float(cov[0]/cov[-1]),'scale':float(np.sqrt((a*a).sum(1).mean()))})
(root/'anderson-v2/metrics.json').write_text(json.dumps(rows,indent=2))
camera=np.array([[1,-1,0],[1,1,-2]],float);camera/=np.linalg.norm(camera,axis=1)[:,None]
for page in range((len(ids)+15)//16):
    fig,axs=plt.subplots(4,4,figsize=(14,14),facecolor='#080e19')
    for ax in axs.flat:ax.set_axis_off()
    for ax,row in zip(axs.flat,range(page*16,min(len(ids),(page+1)*16))):
        i=ids[row];z=tails[row];xy=(z-z.mean(0))@camera.T
        scale=np.max(np.abs(xy))*1.04
        hist,_,_=np.histogram2d(xy[:,0],xy[:,1],bins=360,range=[[-scale,scale]]*2)
        ax.imshow(np.log1p(hist.T),origin='lower',cmap='magma')
        ax.set_title(f'{i}: kind {int(d["parameters"][i,3])}, k={d["parameters"][i,0]:.3f}\nL={d["windows"][i].mean():.3f}',color='white',fontsize=10)
    fig.tight_layout();fig.savefig(root/f'anderson-v2/gallery-{page}.png',dpi=110,facecolor=fig.get_facecolor());plt.close(fig)
print(json.dumps({'positive_count':len(ids),'by_kind':{str(k):int((d['parameters'][ids,3]==k).sum()) for k in [0,2,3,4]}}))

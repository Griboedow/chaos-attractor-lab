"""External diagnostic contact sheet; finite windows, not attractor portraits."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

panels=[]
for n in (3,4):
    data=json.loads(Path(f'followup01/n{n}.json').read_text())
    a=np.load(f'followup01/n{n}-final-25000-steps.npz')
    for r,z in zip(data['runs'],a['states']):
        if r['blocks'][-1]['rate']>.01:panels.append((n,r,z))
fig,axs=plt.subplots(4,4,figsize=(12,12),constrained_layout=True)
for ax,(n,r,z) in zip(axs.flat,panels):
    centered=z-z.mean(0);u,s,vh=np.linalg.svd(centered,full_matrices=False);xy=centered@vh[:2].T
    ax.hexbin(*xy.T,gridsize=120,bins='log',mincnt=1,cmap='Blues',linewidths=0)
    ax.set_title(f"{n}x{n} Case {r['case']}; rate {r['blocks'][-1]['rate']:.3f}",fontsize=10)
    ax.set_aspect('equal');ax.set_xlabel('PC1');ax.set_ylabel('PC2')
for ax in axs.flat[len(panels):]:ax.axis('off')
fig.suptitle('Final 25,000 iterations at step 530,000; unverified finite windows',fontsize=13)
fig.savefig('preview-530k.png',dpi=130)

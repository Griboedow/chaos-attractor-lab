from pathlib import Path
import sys
from icosa_paths import ROOT,RESEARCH,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=OUTPUT;out=root/'icosa-v1';d={**np.load(out/'input.npz'),**np.load(out/'data.npz')};ids=np.array(json.loads((out/'summary.json').read_text())['positive_ids']);rows=[]
rng=np.random.default_rng(2026110145)
for i in ids:
    z=d['tail'][i];rms=np.sqrt((z*z).sum(1).mean());tree=cKDTree(z[:2000]);query=z[2000:][::2]
    baseline=float(np.sqrt(np.mean(tree.query(query)[0]**2))/rms)
    errors=[float(np.sqrt(np.mean(tree.query(query@g.T)[0]**2))/rms) for g in d['group']]
    L=float(d['windows'][i].mean());rows.append({'id':int(i),'L':L,'minL':float(d['windows'][i].min()),'gamma':float(d['gamma'][i]),'sigma':float(d['sigma'][i]),'estimated_KY':float(2+L/(3*d['gamma'][i]+L)),'I_self_ratio':max(errors)/baseline,'I_error':max(errors),'self_error':baseline,'mean_over_rms':float(np.linalg.norm(z.mean(0))/rms)})
rows.sort(key=lambda r:r['estimated_KY']);(out/'metrics.json').write_text(json.dumps(rows,indent=2))
az=.4;el=.6;camera=np.array([[np.cos(az),np.sin(az),0],[-np.sin(az)*np.sin(el),np.cos(az)*np.sin(el),np.cos(el)]])
for page in range((len(rows)+15)//16):
    fig,axs=plt.subplots(4,4,figsize=(14,14),facecolor='#080e19')
    for ax in axs.flat:ax.set_axis_off()
    for ax,row in zip(axs.flat,rows[page*16:page*16+16]):
        xy=d['tail'][row['id']]@camera.T;half=np.abs(xy).max()*1.04;h,*_=np.histogram2d(*xy.T,bins=400,range=[[-half,half]]*2)
        ax.imshow(np.log1p(h.T),origin='lower',cmap='magma');ax.set_title(f"{row['id']}: L={row['L']:.3f}, D~{row['estimated_KY']:.2f}\nI/self={row['I_self_ratio']:.1f}",color='white',fontsize=10)
    fig.tight_layout();fig.savefig(out/f'gallery-{page}.png',dpi=110);plt.close(fig)
print(json.dumps({'most_symmetric':[{'id':r['id'],'L':r['L'],'D_proxy':r['estimated_KY'],'I_self_ratio':r['I_self_ratio']} for r in sorted(rows,key=lambda r:r['I_self_ratio'])[:10]]}),flush=True)

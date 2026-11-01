from paths import OUTPUT, ROOT, RESEARCH
from pathlib import Path
import json
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=OUTPUT
pilot=np.load(root/'tetra-v1/data.npz');d=np.load(root/'tetra-long/data.npz')
az=.4;el=.6;camera=np.array([[np.cos(az),np.sin(az),0],[-np.sin(az)*np.sin(el),np.cos(az)*np.sin(el),np.cos(el)]])
rows=[];rng=np.random.default_rng(2026110135)
for lane,i in enumerate(d['ids']):
    z=d['tail'][lane];rms=np.sqrt((z*z).sum(1).mean());reference=z[:25000][rng.choice(25000,10000,False)];query=z[25000:][rng.choice(25000,3000,False)]
    tree=cKDTree(reference);errors=[float(np.sqrt(np.mean(tree.query(query@g.T)[0]**2))/rms) for g in pilot['group']]
    baseline=float(np.sqrt(np.mean(tree.query(query)[0]**2))/rms)
    rows.append({'id':int(i),'lane':lane,'L':float(d['windows'][lane].mean()),'minL':float(d['windows'][lane].min()),'gamma':float(pilot['gamma'][i]),'T_self_ratio':max(errors)/baseline,'T_errors':errors,'self':baseline,'mean':z.mean(0).tolist()})
(root/'tetra-long/metrics.json').write_text(json.dumps(rows,indent=2))
for page in range((len(rows)+11)//12):
    fig,axs=plt.subplots(3,4,figsize=(16,12),facecolor='#080e19')
    for ax in axs.flat:ax.set_axis_off()
    for ax,r in zip(axs.flat,rows[page*12:page*12+12]):
        z=d['tail'][r['lane']];xy=z@camera.T;scale=np.abs(xy).max()*1.04
        h,*_=np.histogram2d(*xy.T,bins=500,range=[[-scale,scale]]*2)
        ax.imshow(np.log1p(h.T),origin='lower',cmap='magma');ax.set_title(f"{r['id']}: L={r['L']:.4f}, min={r['minL']:.4f}\nT/self={r['T_self_ratio']:.2f}",color='white')
    fig.tight_layout();fig.savefig(root/f'tetra-long/gallery-{page}.png',dpi=110);plt.close(fig)
print(json.dumps([{'id':r['id'],'L':r['L'],'T_self_ratio':r['T_self_ratio']} for r in rows if r['minL']>.003]))

d=np.load(root/'anderson-long/data.npz')
ids=[41,54,73,102,165,449,555,616,629,664,678,682,1178,1218,1538,1616,2076,2179,2268,2332,2524,2716,3036,3043]
for page in range(2):
    fig,axs=plt.subplots(3,4,figsize=(16,12),facecolor='#080e19')
    for ax in axs.flat:ax.set_axis_off()
    for ax,i in zip(axs.flat,ids[page*12:page*12+12]):
        lane=int(np.where(d['ids']==i)[0][0]);z=d['tail'][lane];xy=(z-z.mean(0))@camera.T;scale=np.abs(xy).max()*1.04
        h,*_=np.histogram2d(*xy.T,bins=500,range=[[-scale,scale]]*2)
        ax.imshow(np.log1p(h.T),origin='lower',cmap='magma');ax.set_title(f"{i}: L={d['windows'][lane].mean():.4f}, k={d['parameters'][lane,0]:.3f}",color='white')
    fig.tight_layout();fig.savefig(root/f'anderson-long/gallery-{page}.png',dpi=110);plt.close(fig)

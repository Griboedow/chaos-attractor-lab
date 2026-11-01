from paths import OUTPUT, ROOT, RESEARCH
from pathlib import Path
import json
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=OUTPUT;d=np.load(base/'tetra-v1/data.npz')
ids=np.where((d['windows'].min(1)>.003)&(d['speed']>1e-4))[0]
metrics=[]
for i in ids:
    z=d['tail'][i];rms=np.sqrt(np.mean(np.sum(z*z,axis=1)));tree=cKDTree(z[::2])
    baseline=float(np.sqrt(np.mean(tree.query(z[1::2])[0]**2))/rms)
    group_errors=[]
    for r in d['group']:
        group_errors.append(float(np.sqrt(np.mean(tree.query(z[1::2]@r.T)[0]**2))/rms))
    ev=np.linalg.eigvalsh(np.cov(z.T));mean=float(np.linalg.norm(z.mean(0))/rms)
    metrics.append({'id':int(i),'gamma':float(d['gamma'][i]),'sigma':float(d['sigma'][i]),
                    'L':float(d['windows'][i].mean()),'minL':float(d['windows'][i].min()),
                    'estimated_KY':float(2+d['windows'][i].mean()/(3*d['gamma'][i]+d['windows'][i].mean())),
                    'mean_ratio':mean,'cov_min_max':float(ev[0]/ev[-1]),
                    'symmetry_max':max(group_errors),'self_rms':baseline,
                    'symmetry_ratio':max(group_errors)/baseline})
metrics.sort(key=lambda x:x['estimated_KY'])
(base/'tetra-v1/metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
az=.4;el=.6;camera=np.array([[np.cos(az),np.sin(az),0],[-np.sin(az)*np.sin(el),np.cos(az)*np.sin(el),np.cos(el)]])
for page in range((len(metrics)+15)//16):
    fig,axes=plt.subplots(4,4,figsize=(13,13),facecolor='#060b17')
    for ax in axes.flat:ax.set_axis_off()
    for ax,row in zip(axes.flat,metrics[page*16:page*16+16]):
        z=d['tail'][row['id']];xy=z@camera.T;half=float(np.abs(xy).max()*1.05)
        h,*_=np.histogram2d(*xy.T,bins=420,range=[[-half,half]]*2)
        ax.imshow(np.log1p(h.T),origin='lower',cmap='magma');ax.set_aspect('equal')
        ax.set_title(f"{row['id']} L={row['L']:.3f} D~{row['estimated_KY']:.2f}\nT/self={row['symmetry_ratio']:.1f}",color='white',fontsize=10)
    fig.tight_layout();fig.savefig(base/f'tetra-v1/gallery-{page}.png',dpi=110);plt.close(fig)
print(json.dumps({'most_thin':metrics[:12],'most_symmetric':sorted(metrics,key=lambda r:r['symmetry_ratio'])[:8]},indent=2))

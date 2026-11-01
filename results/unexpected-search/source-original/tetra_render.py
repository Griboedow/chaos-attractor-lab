"""Orthographic portraits from one trajectory; no symmetry copying."""
from pathlib import Path
import json,time,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from tetra_vortex import survey
sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')
from render_twist import density,colorize

root=Path(__file__).parent;d=np.load(root/'tetra-v1/data.npz');out=root/'tetra-render';out.mkdir(exist_ok=True)
ids=np.array([421,565]);begin=time.time()
w,clouds,ends,speed,radius=survey(d['centers'][ids],d['axes'][ids],d['sigma'][ids],d['gamma'][ids],d['starts'][ids],240000,3200000,.0125,3200000)
np.savez_compressed(out/'dynamics.npz',ids=ids,windows=w,end=ends,speed=speed,radius=radius)
for lane,ident in enumerate(ids):np.save(out/f'dense-{ident}.npy',clouds[lane])
rows=[]
angles=[(0,0),(.4,.6),(.785398,.61548),(.9,.15)]
fig,axs=plt.subplots(2,4,figsize=(18,9),facecolor='#080e19')
for lane,ident in enumerate(ids):
    points=clouds[lane]
    np.savez_compressed(out/f'sample-{ident}.npz',points=points[::32])
    for j,(az,el) in enumerate(angles):
        camera=np.array([[np.cos(az),np.sin(az),0],[-np.sin(az)*np.sin(el),np.cos(az)*np.sin(el),np.cos(el)]])
        xy=points@camera.T;half=float(np.abs(xy).max()*1.04);hist,_,_=np.histogram2d(*xy.T,bins=1000,range=[[-half,half]]*2)
        rgb,_=colorize(hist.T[::-1])
        rgb.save(out/f'{ident}-view-{j}.png')
        ax=axs[lane,j];ax.imshow(rgb);ax.set_axis_off();ax.set_title(f'{ident}, camera {j}',color='white')
        rows.append({'id':int(ident),'view':j,'camera':camera.tolist(),'half_width':half,'points':len(points),'warm_time':3000,'observed_time':40000,'h':.0125,'windows':w[lane].tolist(),'speed':float(speed[lane]),'radius':float(radius[lane])})
    print(json.dumps({'id':int(ident),'windows':w[lane].tolist(),'speed':float(speed[lane])}),flush=True)
fig.tight_layout();fig.savefig(out/'views.png',dpi=100);plt.close(fig)
(out/'render.json').write_text(json.dumps({'seconds':time.time()-begin,'rows':rows},indent=2))

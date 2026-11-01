from pathlib import Path
import json,sys
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')
from render_twist import colorize

root=Path(__file__).parent;d=np.load(root/'icosa-v1/input.npz');v=np.load(root/'icosa-spectra/data.npz');rows=[]
fig,axs=plt.subplots(4,4,figsize=(16,16),facecolor='#080e19')
for lane,ident in enumerate([29,183,43,177]):
    idx=int(np.where((v['ids']==ident)&(v['h']==.00625))[0][0]);z=v['tail'][idx];rms=np.sqrt((z*z).sum(1).mean());tree=cKDTree(z[:20000]);query=z[20000:][::10]
    errors=[float(np.sqrt(np.mean(tree.query(query@g.T)[0]**2))/rms) for g in d['group']];identity=float(np.sqrt(np.mean(tree.query(query)[0]**2))/rms)
    row={'id':ident,'max_rotated_error':max(errors),'identity_error':identity,'group_over_identity':max(errors)/identity,'mean_over_rms':float(np.linalg.norm(z.mean(0))/rms),'spectrum':v['windows'][idx].mean(0).tolist()};rows.append(row)
    for j,(az,el) in enumerate([(0,0),(.4,.6),(.8,.6),(.9,.15)]):
        camera=np.array([[np.cos(az),np.sin(az),0],[-np.sin(az)*np.sin(el),np.cos(az)*np.sin(el),np.cos(el)]])
        xy=z@camera.T;half=np.abs(xy).max()*1.04;h,*_=np.histogram2d(*xy.T,bins=700,range=[[-half,half]]*2);pic,_=colorize(h.T[::-1])
        ax=axs[lane,j];ax.imshow(pic);ax.set_axis_off();ax.set_title(f"{ident}, view {j}; I/self={row['group_over_identity']:.2f}",color='white')
fig.tight_layout();fig.savefig(root/'icosa-spectra/gallery.png',dpi=100);plt.close(fig)
(root/'icosa-spectra/symmetry.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows),flush=True)

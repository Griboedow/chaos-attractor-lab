"""Plot actual numerical periodic candidates against the archived chaotic cloud."""
from pathlib import Path
import sys
from cycle_paths import ROOT,RESEARCH,INPUT,STRUCTURE,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(RESEARCH))
from render_twist import density,colorize

base=OUTPUT;out=base/'icosa-cycle-visuals-final';out.mkdir(exist_ok=False)
z=np.load(INPUT/'icosa-long/177-dense.npy',mmap_mode='r');foci=np.load(STRUCTURE/'icosa-itinerary-final/data.npz')['foci']
records=json.loads((base/'icosa-cycle-refinement/report.json').read_text())['rows'];it=json.loads((base/'icosa-directed-cycles/report.json').read_text())['cycles']
fig,axs=plt.subplots(1,2,figsize=(12,6.4),facecolor='#070c15');metadata=[]
for ax,row,seq,color,label in zip(axs,records,it,['#65e5ff','#ffc46d'],['Triangle circuit','Pentagon circuit']):
    points=np.load(base/f"icosa-cycle-refinement/cycle-{row['index']}.npz")['points'];ids=seq['visits_cyclic'];forward=foci[ids].mean(0);forward/=np.linalg.norm(forward);right=foci[ids[0]]-forward*(foci[ids[0]]@forward);right/=np.linalg.norm(right);up=np.cross(forward,right);cam=np.array([right,up,forward]);half=1.36
    bg,_=colorize(density(z,cam[:2],1200,half));ax.imshow((np.asarray(bg)*.42).astype(np.uint8),extent=[-half,half,-half,half]);xy=points@cam[:2].T
    ax.plot(*xy.T,lw=1.15,color=color);mark=foci[ids]@cam[:2].T;ax.scatter(*mark.T,s=23,facecolors='none',edgecolors='white',linewidths=.65)
    for k in [2000,13000,26000]:
        a=xy[k];b=xy[k+120];ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':'white','lw':.9})
    ax.set(xlim=(-half,half),ylim=(-half,half),aspect='equal');ax.set_axis_off();ax.set_title(f"{label} / T = {row['period']:.6f}",color='white',fontsize=13)
    metadata.append({'index':row['index'],'camera':cam.tolist(),'half_width':half,'states':len(points),'arrow_sample_indices':[2000,13000,26000],'marked_focus_ids':ids})
fig.text(.5,.07,'Colored curves: two unstable numerical periodic orbits of case 177',ha='center',color='#dce7fa',fontsize=11)
fig.text(.5,.035,'Dim background: the actual chaotic cloud. Circles and arrows are annotations.',ha='center',color='#9dafc8',fontsize=10)
fig.tight_layout(rect=[0,.09,1,.95]);fig.savefig(out/'icosahedral-cycles.png',dpi=180,bbox_inches='tight',pad_inches=.18);plt.close(fig)
(out/'render.json').write_text(json.dumps({'panels':metadata,'limits':'Separate rigid cameras align each circuit with the viewer. Colored curves are direct full-period integrations, not rotated copies. Periodic existence is not interval certified.'},indent=2));print('Cycle figure complete',flush=True)

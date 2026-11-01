"""Explicitly annotated structural diagrams and actual-time trajectory animation."""
from pathlib import Path
import sys
from structure_paths import ROOT,RESEARCH,INPUT,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image,ImageDraw
sys.path.insert(0,str(RESEARCH))
from render_twist import density,colorize

base=OUTPUT;out=base/'icosa-mechanism-visuals';out.mkdir(exist_ok=False)
z=np.load(INPUT/'icosa-long/177-dense.npy',mmap_mode='r');raw=np.load(INPUT/'icosa-final/portrait-177.npz');cam=raw['camera'];half=float(raw['half_width']);graph=np.load(base/'icosa-itinerary-final/data.npz');foci=graph['foci'];counts=graph['transition_counts'];xy=foci@cam[:2].T
background,_=colorize(density(z,cam[:2],1000,half));fig,axs=plt.subplots(1,2,figsize=(12,6),facecolor='#080e19')
axs[0].imshow(background,extent=[-half,half,-half,half]);axs[0].scatter(*xy.T,s=15,c='#55e8ff',edgecolors='white',linewidths=.3)
segments=[];weights=[]
for i in range(30):
    for j in range(i):
        if counts[i,j]+counts[j,i]>0:segments.append([xy[i],xy[j]]);weights.append(int(counts[i,j]+counts[j,i]))
axs[1].add_collection(LineCollection(segments,colors='#8e62ce',linewidths=.4+2*np.array(weights)/max(weights),alpha=.8));axs[1].scatter(*xy.T,s=25,c='#55e8ff',edgecolors='white',linewidths=.4)
for ax,title in zip(axs,['Actual orbit + 30 saddle-focus centers','Observed transitions, projected onto the same camera']):
    ax.set_facecolor('black');ax.set_xlim(-half,half);ax.set_ylim(-half,half);ax.set_aspect('equal');ax.set_axis_off();ax.set_title(title,color='white',fontsize=11)
fig.text(.5,.03,'Cyan points and straight guide edges are annotations, not added trajectory states.',ha='center',color='#bdc8dc',fontsize=10)
fig.tight_layout(rect=[0,.05,1,.98]);fig.savefig(out/'icosahedral-skeleton.png',dpi=160);plt.close(fig)
resolution=680;bg,_=colorize(density(z[::4],cam[:2],resolution,half),vmax=np.log1p(80));bg=Image.fromarray((np.asarray(bg)*.35).astype(np.uint8))
coords=z[:25600]@cam[:2].T;pixels=np.column_stack([(coords[:,0]/half+1)*resolution/2,(1-coords[:,1]/half)*resolution/2]);frames=[];sample_indices=np.linspace(0,len(coords)-1,160,dtype=int)
for index in sample_indices:
    frame=Image.new('RGB',(resolution,resolution+42),'black');frame.paste(bg,(0,20));draw=ImageDraw.Draw(frame)
    trail=pixels[max(0,index-960):index+1:2]+[0,20]
    if len(trail)>1:draw.line([tuple(p) for p in trail],fill='#bfefff',width=1)
    x,y=pixels[index]+[0,20];draw.ellipse([x-3,y-3,x+3,y+3],fill='#ffffff')
    draw.text((12,4),f'Case 177 / actual time {2000+index*.0125:.2f}',fill='#d9e8fa')
    draw.text((12,resolution+25),'Fixed camera; actual moving state; dim full-cloud context',fill='#91a5bd');frames.append(frame)
frames[0].save(out/'icosahedral-motion.gif',save_all=True,append_images=frames[1:],duration=80,loop=0,optimize=False)
(out/'render.json').write_text(json.dumps({'id':177,'camera':cam.tolist(),'half_width':half,'diagram':'The cyan equilibrium markers and straight edge guides are explicit annotations. Edge width encodes the sum of observed transitions in both directions.','animation':{'frames':160,'resolution':resolution,'initial_ode_time':2000,'time_per_stored_state':.0125,'sample_indices':sample_indices.tolist(),'trail_states':960,'background_states':320000,'background_brightness_factor':.35,'loop_note':'The animation restarts at its beginning; that playback jump is not a dynamical recurrence.'}},indent=2))
print('Structural figure and actual-time animation complete',flush=True)

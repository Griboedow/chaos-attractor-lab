"""Render recorded local returns; connecting lines are only visual guides."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

base=Path(__file__).parent;out=base/'icosa-return-visuals';out.mkdir(exist_ok=False)
d=np.load(base/'icosa-return-grid/data.npz');x=d['inputs'];y=d['outputs'];strips=d['strips'];V=float(d['stable_half_width'])
colors=['#007f90','#b45514']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold'})
fig=plt.figure(figsize=(14,9),facecolor='#fcfbf7');gs=fig.add_gridspec(2,2,width_ratios=[1,1.22],hspace=.48,wspace=.3)
for k in range(2):
    ax=fig.add_subplot(gs[k,0]);ax.set_facecolor('#fcfbf7')
    for lo,hi in strips:ax.axhspan(lo*1000,hi*1000,color='#484e56',alpha=.17,zorder=0)
    for j in range(17):
        ax.plot(x[k,:,j,0]*1e6,y[k,:,j,0]*1000,color=colors[k],alpha=.16,lw=.8)
        ax.scatter(x[k,:,j,0]*1e6,y[k,:,j,0]*1000,color=colors[k],s=5,alpha=.5)
    ax.set(title=f'Source strip {k}: sampled expansion',xlabel='Initial u (millionths)',ylabel='Returned u (thousandths)')
    ax.grid(alpha=.13);ax.text(.02,.02,'Gray bands = the two target u intervals',transform=ax.transAxes,fontsize=9,color='#444444')
ax=fig.add_subplot(gs[:,1]);ax.set_facecolor('#fcfbf7')
for k,(lo,hi) in enumerate(strips):
    ax.add_patch(Rectangle((lo*1000,-V*1e6),(hi-lo)*1000,2*V*1e6,facecolor=colors[k],edgecolor=colors[k],lw=1.5,alpha=.4))
    ax.text((lo+hi)*500,218,f'Target {k}',ha='center',color=colors[k],weight='bold')
for k in range(2):
    for j in range(17):ax.plot(y[k,:,j,0]*1000,y[k,:,j,1]*1e6,color=colors[k],lw=.55,alpha=.24)
    ax.scatter(y[k,:,:,0].ravel()*1000,y[k,:,:,1].ravel()*1e6,s=9,color=colors[k],alpha=.72,label=f'Image of strip {k}')
ax.axhline(200,c='#63676b',ls='--',lw=.8);ax.axhline(-200,c='#63676b',ls='--',lw=.8)
ax.set(xlabel='Returned u (thousandths)',ylabel='Returned v (millionths)',title='Return images inside the stable bounds',ylim=(-245,255))
ax.grid(alpha=.13);ax.legend(loc='lower left',frameon=False)
fig.suptitle('Two candidate return strips inside the icosahedral flow',fontsize=19,weight='bold',y=.97)
fig.text(.5,.922,'Case 177  |  17 x 17 initial points per strip  |  Adaptive ODE and event-time derivative',ha='center',fontsize=11,color='#444444')
fig.text(.5,.035,'Dots are computed returns; lines join samples. No all-point covering, horseshoe or entropy bound is certified.',ha='center',fontsize=10,color='#444444')
fig.subplots_adjust(top=.84,bottom=.12,left=.075,right=.975)
fig.savefig(out/'icosahedral-returns.png',dpi=180,facecolor=fig.get_facecolor());plt.close(fig)
(out/'render.json').write_text(json.dumps({'source':'icosa-return-grid/data.npz','image':'icosahedral-returns.png','sample_count':578,'panels':'Two input-output graphs and one returned-state section, all from the stored grid.','connecting_lines':'Visual interpolation only; no enclosure between samples.','transforms':{'u_input':1e6,'u_output':1000,'v_output':1e6},'limits':'A diagnostic of an existing candidate, not a new attractor portrait or proof of chaos.'},indent=2))

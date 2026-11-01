"""Inspect the scout's crossings under narrower, explicitly declared patches."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base=Path(__file__).parent;out=base/'icosa-return-analysis';out.mkdir(exist_ok=False);d=np.load(base/'icosa-return-scout/data.npz');e=d['events'];a=d['parameters'];rows=[];selected={}
for U,V in [(.03,.003),(.1,.005),(.15,.01)]:
    good=(abs(e[:,:,1])<U)&(abs(e[:,:,2])<V);found=good.any(1);which=good.argmax(1);first=e[np.arange(len(a)),which].copy();first[~found]=np.nan
    rows.append({'u_half_width':U,'v_half_width':V,'returns':int(found.sum()),'times_percentiles':np.nanpercentile(first[:,0],[0,25,50,75,100]).tolist()});selected[f'{U}-{V}']=first
first=selected['0.1-0.005'];crossings=[]
for i in range(len(a)-1):
    if np.isfinite(first[i:i+2]).all() and first[i,1]*first[i+1,1]<0 and abs(first[i+1,0]-first[i,0])<2:
        crossings.append({'indices':[i,i+1],'input':[float(a[i]),float(a[i+1])],'output_u':first[i:i+2,1].tolist(),'output_v':first[i:i+2,2].tolist(),'times':first[i:i+2,0].tolist()})
fig,axs=plt.subplots(1,3,figsize=(12,4));sc=axs[0].scatter(a,first[:,1],c=first[:,0],s=7,cmap='viridis');axs[0].plot(a,a,c='gray',lw=.5);axs[0].axhline(0,c='gray',lw=.5);axs[0].set(xlabel='Initial unstable coordinate',ylabel='First local return: unstable coordinate')
axs[1].scatter(a,first[:,2],c=first[:,0],s=7,cmap='viridis');axs[1].set(xlabel='Initial unstable coordinate',ylabel='Return stable coordinate');axs[2].scatter(first[:,1],first[:,2],c=first[:,0],s=7,cmap='viridis');axs[2].set(xlabel='Return unstable coordinate',ylabel='Return stable coordinate');fig.colorbar(sc,ax=axs[2],label='Return time');fig.suptitle('Finite RK4 scout / local patch |u| < 0.1, |v| < 0.005');fig.tight_layout();fig.savefig(out/'preview.png',dpi=150);plt.close(fig)
report={'patches':rows,'sampled_sign_crossing_candidates':crossings,'limits':'A coarse sampled line, not an invariant manifold or a covering proof. Return branches may be missed or switch between samples.'};(out/'report.json').write_text(json.dumps(report,indent=2));np.savez_compressed(out/'first-returns.npz',parameters=a,**selected);print(json.dumps(report),flush=True)

"""External diagnostic of the three slender-rod scouts; not an attractor portrait."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import spatial_fast as model

root=Path(__file__).resolve().parent
folders=[p for p in sorted(root.iterdir()) if (p/'trajectory.npz').exists()]
fig=plt.figure(figsize=(13,4*len(folders)))
for row,p in enumerate(folders):
    d=np.load(p/'trajectory.npz');t=d['t'];q=d['q'];tip=d['tip'];late=t>=max(t[-1]/2,t[-1]-20)
    r=json.loads((p/'report.json').read_text(encoding='utf-8'))
    ax=fig.add_subplot(len(folders),3,row*3+1,projection='3d')
    for index in np.linspace(np.flatnonzero(late)[0],len(t)-1,12,dtype=int):
        xyz=model.geometry(q[index])[1];ax.plot(*xyz.T,color='#1e6091',alpha=.45,lw=.9)
    ax.set(xlabel='x',ylabel='y',zlabel='z',title=f"Field {r['parameters']['chi']:g}: late shapes")
    ax.set_box_aspect((1,1,1))
    ax=fig.add_subplot(len(folders),3,row*3+2)
    ax.plot(tip[late,0],tip[late,1],color='#1e6091',lw=.6)
    ax.set(xlabel='tip x',ylabel='tip y',title='Late tip trace');ax.set_aspect('equal',adjustable='datalim')
    ax=fig.add_subplot(len(folders),3,row*3+3)
    ax.plot(t,tip[:,2],color='#1e6091',lw=.6)
    ax.axvline(t[late][0],color='0.5',ls='--',lw=.6)
    ax.set(xlabel='time',ylabel='tip height',title='Completed' if r['success'] else 'Newton failure')
fig.tight_layout();fig.savefig(root/'diagnostic.png',dpi=150);plt.close(fig)

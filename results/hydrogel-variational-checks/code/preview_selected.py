"""External diagnostics; no attractor or novelty claim."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import spatial_accelerated as model

root=Path(__file__).resolve().parent
names=['spectrum117-n64','spectrum124-n64','spectrum127-n64','field130-continuation']
fig=plt.figure(figsize=(12,10))
for row,name in enumerate(names):
    with np.load(root/name/'trajectory.npz') as d:
        t=d['t'];q=d['q'];tip=d['tip']
    report=json.loads((root/name/'report.json').read_text(encoding='utf-8'))
    late=t>=t[-1]-20
    ax=fig.add_subplot(len(names),3,row*3+1,projection='3d')
    for i in np.linspace(np.flatnonzero(late)[0],len(t)-1,10,dtype=int):
        xyz=model.geometry(q[i])[1];ax.plot(*xyz.T,color='#1e6091',alpha=.5,lw=.8)
    ax.set(xlabel='x',ylabel='y',zlabel='z',title=f"Field {report['parameters']['chi']:g}: late shapes")
    ax.set_box_aspect((1,1,1))
    ax=fig.add_subplot(len(names),3,row*3+2)
    ax.plot(tip[late,0],tip[late,1],color='#1e6091',lw=.5)
    ax.set(xlabel='Tip x',ylabel='Tip y',title='Last 20 time units');ax.set_aspect('equal',adjustable='datalim')
    ax=fig.add_subplot(len(names),3,row*3+3)
    ax.plot(t,tip[:,2],color='#1e6091',lw=.7)
    ax.set(xlabel='Time',ylabel='Tip height',title='Complete' if report['success'] else 'Newton failure')
fig.tight_layout();fig.savefig(root/'diagnostic-selected.png',dpi=130);plt.close(fig)

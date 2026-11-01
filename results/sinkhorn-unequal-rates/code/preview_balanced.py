"""Rank finite-window spatial variance; render orthogonal projections without whitening."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
data=json.loads((args.input/'3x2.json').read_text());arr=np.load(args.input/'3x2-final-50000-states.npz');rows=[]
for r,z in zip(data['runs'],arr['states']):
    if not np.isfinite(z).all():continue
    mean=z.mean(0);u,s,vh=np.linalg.svd(z-mean,full_matrices=False);variance=s*s/np.sum(s*s)
    rows.append(dict(case=r['case'],growth=r['blocks'][-1]['rate'],variance=variance.tolist(),rotation=vh.tolist(),center=mean.tolist(),states=z))
rows.sort(key=lambda r:r['variance'][2],reverse=True)
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9});fig=plt.figure(figsize=(14,10),layout='constrained')
for i,r in enumerate(rows[:12]):
    q=(r['states'][::5]-r['center'])@np.array(r['rotation']).T;ax=fig.add_subplot(3,4,i+1,projection='3d');ax.scatter(*q.T,s=.25,c='#164e63',alpha=.3)
    center=(q.max(0)+q.min(0))/2;span=np.max(np.ptp(q,axis=0))/2
    for setlim,c in zip((ax.set_xlim,ax.set_ylim,ax.set_zlim),center):setlim(c-span,c+span)
    ax.set_box_aspect((1,1,1));ax.view_init(22,-55);ax.set_title(f"Case {r['case']}; growth {r['growth']:.3f}\nPC3 variance {r['variance'][2]:.4f}");ax.set_xlabel('PC1');ax.set_ylabel('PC2');ax.set_zlabel('PC3')
fig.suptitle('Balanced unequal-rate kernels: finite windows, attraction unverified');fig.savefig(args.output,dpi=150)
for r in rows:del r['states']
args.output.with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows[:12]))

"""External geometric diagnostic of unequal-rate finite trajectories."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9});fig=plt.figure(figsize=(13,7.5),layout='constrained');panel=1;records=[]
for m,n in [(2,3),(3,2)]:
    data=json.loads((args.input/f'{m}x{n}.json').read_text());a=np.load(args.input/f'{m}x{n}-final-50000-states.npz')
    for r,z in zip(data['runs'],a['states']):
        if not np.isfinite(z).all():continue
        z=z[::5];mean=z.mean(0);u,s,vh=np.linalg.svd(z-mean,full_matrices=False);q=(z-mean)@vh.T
        ax=fig.add_subplot(2,4,panel,projection='3d');panel+=1
        ax.scatter(*q.T,s=.2,c='#164e63',alpha=.25,rasterized=True);ax.view_init(24,-52)
        span=np.max(np.ptp(q,axis=0))/2;center=(q.max(0)+q.min(0))/2
        ax.set_xlim(center[0]-span,center[0]+span);ax.set_ylim(center[1]-span,center[1]+span);ax.set_zlim(center[2]-span,center[2]+span);ax.set_box_aspect((1,1,1))
        ax.set_title(f"{m} × {n}, Case {r['case']}\nrate {r['blocks'][-1]['rate']:.4f}")
        ax.set_xlabel('PC1');ax.set_ylabel('PC2');ax.set_zlabel('PC3')
        records.append(dict(shape=[m,n],case=r['case'],pca_variance=(s*s/np.sum(s*s)).tolist(),rotation=vh.tolist(),center=mean.tolist()))
fig.suptitle('Finite trajectories of unequal-rate Sinkhorn; asymptotic attraction unverified')
fig.savefig(args.output,dpi=150);args.output.with_suffix('.json').write_text(json.dumps(records,indent=2)+'\n')

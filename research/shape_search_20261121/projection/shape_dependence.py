"""Compare joint spatial occupancy with shifted PCA marginals of the same orbit."""
import argparse
import json
from pathlib import Path
import numpy as np


def entropy(keys,spaces):
    n=np.bincount(keys,minlength=spaces)
    p=n[n>0]/len(keys)
    return -np.dot(p,np.log(p))


def multiinfo(labels,k):
    a,b,c=labels.T
    return sum(entropy(v,k) for v in (a,b,c))-entropy((a*k+b)*k+c,k**3)


def measure(z,seed):
    z=z.astype(float)
    z-=z.mean(0)
    vals,vec=np.linalg.eigh(z.T@z)
    coord=z@vec[:,::-1]
    rng=np.random.default_rng(seed)
    rows=[]
    for k in (8,16,24):
        bins=np.quantile(coord,np.arange(1,k)/k,axis=0)
        index=np.column_stack([np.searchsorted(bins[:,j],coord[:,j]) for j in range(3)])
        for half,p in enumerate(np.array_split(index,2)):
            obs=multiinfo(p,k);null=[]
            for _ in range(16):
                offsets=rng.integers(700,len(p)-700,size=2)
                shifted=np.column_stack((p[:,0],np.roll(p[:,1],offsets[0]),np.roll(p[:,2],offsets[1])))
                null.append(multiinfo(shifted,k))
            rows.append(dict(bins=k,half=half+1,observed=float(obs),shifted_median=float(np.median(null)),shifted_max=float(np.max(null)),excess_above_shifted_max=float(obs-max(null))))
    return rows


def main(args):
    rows=[]
    for case in json.loads((args.source/'manifest.json').read_text()):
        with np.load(case['cloud']) as f:cloud=f['points']
        detail=measure(cloud,47213+case['case'])
        rows.append(dict(case=case['case'],metrics=detail))
    doc=dict(method='PCA in raw coordinate units; quantile bins learned on the full orbit; sum of marginal entropies minus joint entropy within each temporal half. Each of 16 controls keeps axis one fixed and independently wraps axes two and three by a random shift of at least 700 points. All entropies use natural logs.',rows=rows)
    args.out.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')
    print('positive all scales and halves',sum(min(x['excess_above_shifted_max'] for x in row['metrics'])>0 for row in rows),'/',len(rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())

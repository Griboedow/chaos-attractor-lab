"""Finite sampled centerline clearance, with an explicit material-gap cutoff."""
import json
from pathlib import Path
import numpy as np
from numba import njit
import spatial_fast as m

@njit(cache=True)
def distance(p,q,a,b):
    w=p-q;aa=np.dot(a,a);ab=np.dot(a,b);bb=np.dot(b,b);aw=np.dot(a,w);bw=np.dot(b,w)
    det=aa*bb-ab*ab
    best=1e300
    if det>1e-16*aa*bb:
        u=(ab*bw-bb*aw)/det;v=(aa*bw-ab*aw)/det
        if 0<=u<=1 and 0<=v<=1:
            d=w+u*a-v*b;best=np.dot(d,d)
    for edge in range(4):
        if edge<2:
            u=float(edge);v=min(1.0,max(0.0,(bw+u*ab)/bb))
        else:
            v=float(edge-2);u=min(1.0,max(0.0,(v*ab-aw)/aa))
        d=w+u*a-v*b;best=min(best,np.dot(d,d))
    return np.sqrt(best)

@njit(cache=True)
def measure(q,gap):
    frames,nodes=m.geometry(q);n=len(q);best=1e300;pair=np.zeros(2,dtype=np.int64)
    for i in range(n):
        for j in range(i+gap,n):
            d=distance(nodes[i],nodes[j],nodes[i+1]-nodes[i],nodes[j+1]-nodes[j])
            if d<best:best=d;pair[0]=i;pair[1]=j
    return best,pair

if __name__=='__main__':
    root=Path(__file__).resolve().parent;rows=[]
    for folder in sorted(root.iterdir()):
        if not folder.is_dir() or not (folder/'report.json').exists() or not (folder/'trajectory.npz').exists():continue
        r=json.loads((folder/'report.json').read_text(encoding='utf-8'))
        d=np.load(folder/'trajectory.npz');t=d['t'];q=d['q'];n=len(q[0]);lam=r['parameters']['slenderness'];eta=r['parameters']['eta']
        gap=max(2,int(np.ceil(4*n/lam))+1);threshold=2*min(1.0,eta)/lam
        stride=max(1,round(0.05/(t[1]-t[0])))
        ids=np.arange(np.searchsorted(t,max(t[-1]/2,t[-1]-20)),len(t),stride)
        values=[];pairs=[]
        for i in ids:
            value,pair=measure(q[i],gap);values.append(value);pairs.append(pair.tolist())
        lowest=int(np.argmin(values))
        row=dict(run=folder.name,samples=len(ids),minimum_segment_index_gap=gap,minimum_material_gap=(gap-1)/n,
                 tube_diameter_threshold=threshold,minimum_clearance=float(values[lowest]),time_of_minimum=float(t[ids[lowest]]),
                 pair_at_minimum=pairs[lowest],samples_below_threshold=int(np.sum(np.array(values)<threshold)),
                 scope='Sampled polygonal centerlines only. The material-gap cutoff omits nearby segments; no between-sample, finite-cross-section or continuum contact certificate.')
        rows.append(row);print(json.dumps(row))
    (root/'clearance.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')

"""Measure full-state spatial diversity and select a reproducible visual inspection set."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist


def descriptors(z):
    z=z.astype(float)
    z-=z.mean(0)
    eig,vec=np.linalg.eigh(z.T@z/len(z))
    if eig[0]<=0 or eig.sum()<1e-7:
        return np.full(42,np.nan)
    z=z@vec[:,::-1]/np.sqrt(eig.sum())
    eig=eig[::-1]/eig.sum()
    rad=np.linalg.norm(z,axis=1)
    pairs=np.linalg.norm(z-np.roll(z,173,axis=0),axis=1)
    # Rotation and uniform-scale invariant (up to PCA degeneracy).
    desc=list(eig[:2])+list(np.quantile(rad,[.05,.20,.40,.65,.85,.98]))+list(np.quantile(pairs,[.05,.20,.40,.65,.85,.98]))
    for i in range(3):
        h=np.sort(z[:,i])
        span=max(h[-1]-h[0],1e-8)
        desc.extend(np.sort(np.diff(h))[::-1][:3]/span)
        desc.extend(np.quantile(np.abs(h),[.25,.75]))
    frac=np.clip((z+2.5)/5,0,1-1e-7)
    for n in (8,16,32):
        u=(frac*n).astype(int)
        keys=(u[:,0]*n+u[:,1])*n+u[:,2]
        _,freq=np.unique(keys,return_counts=True)
        p=freq/len(z)
        desc.extend([np.log(len(freq))/np.log(len(z)),float(-(p*np.log(p)).sum()/np.log(len(z)))])
    ids=np.arange(0,len(z),max(1,len(z)//48))[:48]
    _,near=cKDTree(z).query(z[ids],k=15)
    loc=z[near]
    loc-=loc.mean(axis=1,keepdims=True)
    evals=np.linalg.eigvalsh(np.swapaxes(loc,1,2)@loc)
    frac_local=evals/np.maximum(evals.sum(axis=1,keepdims=True),1e-14)
    desc.extend(np.quantile(frac_local[:,0],[.25,.5,.75]))
    desc.extend(np.quantile(frac_local[:,1],[.25,.5,.75]))
    desc.append(float(np.minimum(rad.max(),8)))
    return np.asarray(desc)


def main(args):
    records=[];desc=[];items=[]
    for path in sorted(args.deep.glob('deep-*.npz')):
        with np.load(path) as f:
            ids=f['ids']; data=f['data']; m=f['metrics'];kind=f['kinds']; points=f['points']
            for j,i in enumerate(ids):
                d=descriptors(points[j])
                if np.isfinite(d).all():
                    desc.append(d)
                    items.append(data[:,j,:].copy())
                    records.append([int(i),int(kind[j]),float(m[:4,j].mean()),float(m[:4,j].min())])
    d=np.asarray(desc); records=np.asarray(records); items=np.asarray(items)
    # Prefer genuinely 3D stable occupancy, with a range of continuous sheets,
    # filaments and separated components. Reject concentrated few-point tails.
    mask=(d[:,0]<.92)&(d[:,1]>.07)&(1-d[:,0]-d[:,1]>.012)&(d[:,33]>.53)&(d[:,-1]<5.5)
    d=d[mask]; records=records[mask]; items=items[mask]
    scaled=(d-np.median(d,axis=0))/np.maximum(np.quantile(d,.90,axis=0)-np.quantile(d,.10,axis=0),.04)
    scaled=np.clip(scaled,-2.5,2.5)
    scaled[:,14:29]*=.70
    scaled[:,29:35]*=.80
    scaled[:,-1]*=.45
    quality=(1-d[:,0])+d[:,33]+d[:,39]+np.minimum(records[:,3],1)*.4-np.maximum(0,d[:,-1]-3)*.20
    selected=[]
    # Seed the gallery with each site geometry. Complete it by farthest-first
    # sampling after trimming extreme outlier descriptors to the above range.
    for k in range(5):
        group=np.where(records[:,1]==k)[0]
        if len(group): selected.append(int(group[np.argmax(quality[group])]))
    best=np.min(cdist(scaled,scaled[selected]),axis=1)
    weights=np.clip(.85+.28*(quality-np.median(quality)),.65,1.15)
    while len(selected)<min(args.limit,len(d)):
        best[selected]=-1
        j=int(np.argmax(best*weights))
        selected.append(j)
        best=np.minimum(best,np.linalg.norm(scaled-scaled[j],axis=1))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    np.savez(args.out,ids=records[selected,0].astype(int),kinds=records[selected,1].astype(int),rates=records[selected,2:],data=items[selected],desc=d[selected],total_eligible=np.array([len(d)]))
    print(json.dumps(dict(screened=len(mask),eligible=int(mask.sum()),selected=len(selected))))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--deep',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--limit',type=int,default=160)
    main(p.parse_args())
